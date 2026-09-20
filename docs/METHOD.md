# How the read-structure call works

## 1. The problem

A ribo-seq read is a footprint wrapped in library scaffolding. Reading 5'→3',
any of these may be present or absent:

```
[5' UMI][RT nt][============ footprint ============][3' UMI][barcode][adapter]
   random  A/T          genomic, variable-length        random   fixed   fixed
```

plus poly(A) tails and template-switch (G/C-rich) motifs. We are given only an
accession or a FASTQ and must report the layout — lengths and identities — with
**no prior knowledge of the protocol**, and refuse when the reads cannot support
an answer.

The output is not the *names* of the parts but what a processing pipeline needs:
how many nt to trim from each end, which stretches are random (so must be used
for UMI deduplication), and which adapter to remove.

## 2. The one idea everything rests on

We do **not** trim before aligning. We align the raw read with
`STAR --alignEndsType Local`, so everything non-genomic — 5' UMI, RT additions,
3' UMI, barcode, adapter — is pushed into the **soft clips** instead of
preventing the alignment. The architecture is then read back out of the clips.

The catch: STAR's clip boundary is **not** where the footprint starts. If the
last base of a 5' UMI happens to match the genome base next to the footprint
(probability 1/4), STAR extends the alignment by one and the clip shrinks.
**Clip lengths leak.**

What does *not* leak is the genomic coordinate the alignment implies for **read
position 0**. Extending the alignment leftwards by one decrements *both* `pos`
and `clip5`, so `pos - clip5` is invariant. Anchor on that, and every read
position `p` maps to a fixed genome base. Then ask the only question that matters:

> Does read base `p` match the genome base implied for it?

* Footprint bases are genomic → they match ~95–99 % of the time.
* UMI / adapter / barcode bases are not genomic → they match ~25 % by chance.

**The transition between those two regimes is the architecture.** Everything
downstream is a way of locating that transition precisely and naming what sits
on either side.

## 4. The quantities (`fqdissect/profile.py`)

Everything is descriptive here — no thresholds, no calling. For each read the
profiler keeps the base string, a per-position **match mask** against the implied
genome base, and the two soft-clip boundaries, then aggregates:

**5'-anchored** (read position 0 = first base; exact, no leakage):
* `p5_match[p]` — fraction of reads whose base `p` matches the genome. This is
  essentially the **CDF of the 5'-clip length**, lifted by a 1/4 chance floor:
  `P(p inside aligned block) + P(p clipped)·¼`.
* `p5_comp[p]` — the A/C/G/T composition at position `p`.
* `clip5_hist` — the 5'-clip-length histogram, i.e. the **survival function**
  that `p5_match` integrates.

**3'-anchored.** Read position is meaningless on the 3' side (footprint length
varies), and the alignment *end* is worse than useless — local alignment stops
at a mismatch, so the base just past it disagrees with the genome *by
construction* (match rate 0.00). The two honest anchors are:
* the **adapter start** — `adap_match` / `adap_comp`, indexed so `A0` is the
  first base of the fixed 3' scaffold and `A-1, A-2, …` walk left toward the
  footprint. The scaffold is found **de novo** (next section).
* the **read's own 3' end** — `t3_match` / `t3_comp`, the fallback when no
  adapter is present (already-trimmed deposits).

**Shape histograms:**
* `footprint_len_hist` — genuine footprints spread across ~26–34 nt; a single
  sharp spike means the "footprints" are a fixed-length artefact (adapter dimer).
* `fpend_to_adapter_gap_hist` — sharp mode = fixed-length 3' construct; broad =
  a variable poly(A) tail sits in between.
* `top5p_locus_frac` — the fraction of reads at the single most common 5'
  genomic coordinate (leakage-invariant). Genuine footprints spread across
  thousands of loci (< ~2 % at any one); a majority at one coordinate is a
  single over-represented species, not footprints.

Reads with indels or splice junctions break the position→genome map and are
skipped (a few %); adapter-dimers and poly(A) reads map to genomic homopolymers
and are dropped by a low-complexity filter on the aligned core.


## 4b. Finding the adapter de novo (`fqdissect/profile.py`)

No adapter sequence is known to fqdissect. Blind adapter detectors — fastp's
single-end mode, DNApi, atropos `detect`, minion — count k-mers in read tails and
assemble the over-represented ones, and they have to guess where the insert ends.
Here the alignment has already said where it ends (`fp_end`: STAR's aligned block,
extended right through chance genome matches), so the search is confined to the
**non-genomic 3' tail** of each read:

1. **Seed.** Count every 12-mer of every tail, once per read, remembering its mean
   distance from the insert end. Drop simple k-mers (near-homopolymers, dinucleotide
   repeats, anything containing a run of ≥ 6 — those are tails and dark cycles). All
   k-mers with at least half the top count then compete on **position**: the one
   *nearest the insert* wins. Frequency alone is the wrong criterion — long reads run
   through a ligated linker straight into the sequencing-primer site and carry both in
   (nearly) every read; the 3' adapter is the one ligated to the insert.
2. **Gate.** The seed must occur in ≥ 10 % of the reads — or in ≥ 300 reads at a *fixed*
   distance from the insert end (gap concentration ≥ 0.30). The second route is the
   long-insert library whose adapter falls off the end of most reads: the minority that
   reaches it still reads the construct out exactly.
3. **Assemble.** Anchor the reads on the seed and take the per-position consensus,
   extending left and right for as long as positions are constant (entropy ≤ 0.9 bit)
   and enough reads cover them. The left walk stops by itself at a UMI (random) or at
   the footprint (variable). A homopolymer run in the left extension is a poly(A) tail
   sitting between footprint and adapter: the scaffold starts after it and sheds leading
   bases of the tail's kind.
4. **Re-anchor** on the assembled scaffold. A match of ≥ 12 nt may begin *before* the
   computed insert end: when the genome happens to continue like the adapter, `fp_end`
   is extended into the adapter and the true start lies inside the "insert".

The result is the fixed 3' **scaffold**: the adapter plus any sample barcode in front
of it (a barcode is constant within one FASTQ, so the data cannot separate the two; the
optional name table can — see the README).

**Abundant species are capped first.** At most `max(10, 0.2 %)` of the loaded reads may
share one 5' locus. One rRNA fragment sequenced ten thousand times otherwise turns its
own last bases into a "constant block" in front of the adapter and its genome context
into the match profile. The structure is a property of every molecule alike, so nothing
is lost, and no contaminant reference is needed.

**A UMI behind the adapter.** Some kits (QIAseq miRNA) read
`[insert][adapter][UMI 12][RT-primer site]`. Anchored on the adapter this is a constant
block, a run of 6–16 uniformly random positions, and a *second* constant block (required
— random-looking positions with nothing constant behind them are just whatever the
reads run into). It is reported as a second UMI, and trimmed outside-in: remove the
second constant block, take the UMI off the exposed 3' end, remove the adapter.

## 5. Reading the 5' side (`call_5prime`)

Two derived quantities do the work, both relative to the sample's own **genomic
plateau** (the median deep-position match rate — a noisy library plateaus at
0.79 where a clean one reaches 0.97, so all thresholds are scaled to it, never
absolute):

* **penetrance**`[p] = (plateau − p5_match[p]) / (plateau − ¼)`, clipped to
  [0,1]. This is the survival function of the per-read 5'-construct length: a
  fixed *L*-nt element sits at penetrance ≈1 for *L* positions then steps to 0;
  an element present in only a fraction *f* of molecules sits at *f*.
* **composition** — entropy and A/T vs G/C bias per position.

The rules, footprint-first:

| element | signature | fate |
|---|---|---|
| **footprint start** | first position that is, and stays, genomic (penetrance < `FOOTPRINT_PEN`) | the boundary; everything left is trimmed |
| **RT untemplated nt** | 1–2 nt against the footprint, at *partial* penetrance **or** A/T-biased | **kept** (it is enzymatic, part of the molecule) |
| **template-switch / linker** | full penetrance, **G/C-biased** (≥ `GC_TEMPLATE`) | trim |
| **UMI** | full penetrance, balanced, high entropy | trim **+ dedup** |
| **barcode** | full penetrance, **fixed base** (low entropy) | trim |

The key discriminations: a **step** in `p5_match` (not a level) marks the UMI
length — a smooth ramp is just a tail of mis-mapped reads and must not be read as
a long UMI. A UMI base is uniform; an **RT base is A/T-biased and partial**,
which is how a non-templated addition is told apart from a UMI base even when the
two abut. **No protocol uses a 1-nt UMI**, so a lone non-genomic 5' base is an
RT addition, not a UMI.

## 6. Reading the 3' side (`call_3prime`)

Walk **left from the adapter start** through `adap_match`, nearest-first:

* Stop at the footprint — a base that is genomic (match ≥ `thr_genomic`) **and
  variable** (entropy > const). The "and variable" matters: a *constant* linker
  base can match the genome well above chance simply because one fixed base is
  compared against a biased genomic neighbourhood, so a high match rate alone
  must not end the walk.
* A real UMI/barcode base is *clearly* non-genomic (match ≈ chance). A position
  with only an *intermediate* match is the footprint's own 3' edge, depressed by
  alignment-edge leakage — not a construct base. If **nothing** between footprint
  and adapter is clearly non-genomic, the walk only caught a footprint edge and
  the 3' construct is dropped (no fabricated 1-nt UMI).
* A **homopolymer** run (poly-A) is constant in composition but variable in
  length, so it is a tail, not a fixed barcode, and is peeled off separately.

What remains between footprint and adapter is segmented by entropy into
`random` (UMI), `degenerate`, and `const` (barcode) blocks, in insert-first
order.

## 7. The four functional categories

The named fields (`umi5`, `barcode3`, …) are re-cast into what actually drives
processing — **origin and fate** — because the residual "±1 UMI-vs-spacer"
ambiguity lives *inside* one category and therefore changes none of the derived
quantities:

| category | penetrance | composition | genome match | fate |
|---|---|---|---|---|
| footprint | ~0 | genomic | high | **keep** — the data |
| random-templated (UMI / spacer) | ~1 | balanced | chance | trim **+ dedup** |
| fixed-templated (barcode, linker, adapter) | ~1 | low entropy | chance | trim |
| enzymatic (RT nt / template-switch) | partial (RT) or ~1 (TS) | A/T (RT), G/C (TS) | chance | RT nt **kept**, TS trimmed |

`infer()` emits a `functional` block with the 5'→3' segment list and the two
load-bearing numbers: **`trim_5p`** (nt to clip from the 5' end to reach the
footprint) and **`dedup_umi_len`** (total random-templated content, the mask a
deduplicator needs).

## 8. Refusing to answer

`undetermined` is returned — never a guess — when the profile carries no
readable boundary: too few usable alignments (< 5 000), no genomic plateau
(deep match rate < 0.55, i.e. alignments too noisy to read), a footprint-length
spike above the artefact gate (adapter dimers / fixed contaminants), an adapter
consensus too short to trim on, or a 3' tail that
cannot be segmented. The reason is written into the JSON and the figure title, and
`fqdissect run` exits with status 3. A poly(A) library whose real adapter lies beyond
the tail is *not* refused: the tail is reported and trimmed, the adapter is reported as
not visible.
