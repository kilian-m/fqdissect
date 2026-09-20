"""Turn locally-aligned reads into the positional statistics that reveal the
read structure.

The one idea this module rests on
--------------------------------
For a read aligned in local mode, STAR's soft-clip boundary is *not* a reliable
marker of where the insert starts. If the last base of a 5' UMI happens to match
the genome base next to the insert (probability 1/4), STAR extends the alignment
by one and the clip shrinks. So clip lengths leak.

What does *not* leak is the genomic coordinate that the alignment implies for
read position 0. Extending the alignment leftwards by one decrements both `pos`
and `clip5`, so `pos - clip5` is invariant. Anchoring on that, every read
position p maps to a fixed genomic base, and we can ask the only question that
matters:

    does read base p match the genome base the alignment implies for it?

Insert bases match ~99% of the time. UMI, adapter and barcode bases are not
genomic, so they match ~25% of the time by chance. The transition between those
two regimes is the read structure.

Anchors
-------
5' side: read position 0. The 5' construct has the same length in every read,
so this anchor is exact.

3' side: read position 0 is useless because the insert length varies. The
alignment end is no good either -- STAR extends until it hits a mismatch, so the
base just past it mismatches *by construction*. The two honest anchors are the
**adapter start** and the **read's own 3' end**, and this module profiles around
both.

The adapter is always found de novo
-----------------------------------
No adapter sequence is known to this module. Blind adapter detectors (fastp's
overlap-free mode, DNApi, atropos `detect`, minion) count k-mers in read tails and
assemble the over-represented ones, and have to guess where the insert ends. Here
the alignment has already said where it ends, so the search is restricted to the
**non-genomic 3' tail** of each read:

1. *seed* -- among the frequent complex k-mers in the tails, the one lying nearest
   the insert (`discover_seed`);
2. *anchor* -- locate the seed in every read (partial at the read end allowed);
3. *assemble* -- the per-position consensus around the anchor, extended left and
   right for as long as the position is constant (`assemble_scaffold`): this is the
   fixed 3' scaffold, i.e. the adapter plus any sample barcode in front of it (a
   barcode is constant within one FASTQ, so it is part of the scaffold);
4. *re-anchor* on the assembled scaffold, so the profile is indexed from its first
   base.

A homopolymer (poly-A tail, poly-G dark cycles) is never an adapter: it is counted
separately and only becomes the anchor when nothing else is there.

Everything here is descriptive. Thresholds and calling live in `fqdissect.infer`.
"""
from __future__ import annotations

import os
from collections import Counter
from dataclasses import dataclass

import numpy as np
import pysam

from .utils import LOG

MAXP5 = 24            # read positions profiled from the 5' end
MAXW = 50             # window profiled around the 3' anchors
SEED_K = 12           # k-mer length used to discover the adapter
SEED_NEAR_FRAC = 0.5  # k-mers at least this frequent (vs. the top one) compete on position
SEED_MIN_FRAC = 0.10  # the seed k-mer must occur in this share of the reads ...
SPECIFIC_OVERLAP = 12 # a scaffold match this long may start inside the "insert" (see anchor_on)
LOCUS_CAP_FRAC = 0.002  # no single 5' locus may contribute more than this share of reads
LOCUS_CAP_MIN = 10
MIN_CORE_ENT = 1.2    # bits; aligned cores below this are homopolymer junk
MAX_CORE_BASE = 0.70  # ... as is any core dominated by a single base
SCAFFOLD_ENT = 0.90   # bits; a scaffold position is constant (= infer's ent_const)
HOMOPOLYMER_MIN = 4   # a run this long inside the left extension is a tail, not scaffold
HOMOPOLYMER_PROBE = 20
BASES = "ACGT"

# The seed must be visible in this share of the reads before the 3' profile is
# anchored on it. The caller applies the SAME gate (`min_anchor_frac`) and passes
# its value in, so the two cannot drift apart.
MIN_ANCHOR_FRAC = 0.15

# ... unless the adapter is hiding behind a long insert. When the molecule is about
# as long as the read, only the short-insert minority sequences far enough to reach
# the adapter -- but those reads read the construct out perfectly well. An adapter
# seen in this many reads, at a FIXED distance from the insert end, is real no
# matter how small its share.
MIN_ANCHOR_READS = 300
MIN_GAP_CONC = 0.30

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# base -> 0..3 for the composition counts; anything else (N) -> 4, counted in the
# denominator but in no base's numerator
_CODE = np.full(256, 4, dtype=np.uint8)
_CODE[np.frombuffer(BASES.encode("ascii"), dtype=np.uint8)] = np.arange(4, dtype=np.uint8)


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def max_mm(overlap: int) -> int:
    """cutadapt-style error budget: 0 errors below 9 nt, then ~12%."""
    return int(0.12 * overlap)


def find_adapter(read: str, adapter: str, min_start: int, min_overlap: int = 7) -> int:
    """Leftmost i >= min_start where read[i:] is a prefix of `adapter` (or vice versa).

    Handles the common case where the adapter runs off the end of the read, so
    only its first few bases are visible.
    """
    n, m = len(read), len(adapter)
    for i in range(max(0, min_start), n - min_overlap + 1):
        ov = min(m, n - i)
        if ov < min_overlap:
            break
        budget = max_mm(ov)
        mm = 0
        ok = True
        for a, b in zip(read[i: i + ov], adapter[:ov]):
            if a != b:
                mm += 1
                if mm > budget:
                    ok = False
                    break
        if ok:
            return i
    return -1


def cigar_is_simple(ct) -> bool:
    """True when the read has no I/D/N, so read pos -> genome pos is a shift."""
    return all(op in (0, 4, 7, 8) for op, _ in ct)


def low_complexity(s: str) -> bool:
    """Adapter dimers and poly-A tails map to genomic homopolymers. Their 'insert'
    is meaningless, and they drag every profile toward noise."""
    if not s:
        return True
    counts = [s.count(b) for b in BASES]
    n = len(s)
    if max(counts) / n > MAX_CORE_BASE:
        return True
    p = np.array([c / n for c in counts if c])
    return float(-(p * np.log2(p)).sum()) < MIN_CORE_ENT


def entropy(row) -> float:
    r = np.asarray(row, dtype=float)
    if np.isnan(r).any() or r.sum() <= 0:
        return float("nan")
    r = r[r > 0]
    return float(-(r * np.log2(r)).sum())


@dataclass(slots=True)
class _Read:
    """One usable alignment, reduced to what the profiles need."""

    seq: str             # read as sequenced (5'->3' of the original read)
    n: int               # len(seq)
    codes: np.ndarray    # uint8, 0..3 = ACGT, 4 = other
    matches: np.ndarray  # bool, read base p == implied genome base
    valid: np.ndarray    # bool, position p has an implied genome base at all
    aln_s: int           # first aligned read position (5' clip length)
    aln_e: int           # first read position past STAR's aligned block
    fp_end: int          # ... extended right through chance matches: insert end


def load_reads(bam_path: str, fasta_path: str, max_reads: int):
    """Read the BAM once; keep the read, the match mask and the two boundaries.

    Also count, per used read, the leakage-invariant genomic coordinate of the 5'
    end (`pos - clip5`). Genuine inserts spread across thousands of loci; a large
    share at one coordinate means the reads are a single over-represented species
    (an adapter dimer aligning to an adapter-like locus, a spike-in, a contaminant).
    """
    fa = pysam.FastaFile(fasta_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    recs: list[_Read] = []
    keys: list[tuple] = []
    loc_hist: Counter = Counter()
    n_seen = n_skip_cigar = n_skip_ctx = n_skip_lowcomp = 0

    for aln in bam:
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        n_seen += 1
        if len(recs) >= max_reads:
            break
        ct = aln.cigartuples
        if not ct or not cigar_is_simple(ct):
            n_skip_cigar += 1
            continue

        seq = aln.query_sequence.upper()
        L = len(seq)
        clipL = ct[0][1] if ct[0][0] == 4 else 0
        clipR = ct[-1][1] if ct[-1][0] == 4 else 0

        # every read position must have an implied genome base, so the context
        # has to reach `clip` bases beyond the aligned block on either side
        pad = min(250, L + 10)
        lo, hi = aln.reference_start - pad, aln.reference_end + pad
        if lo < 0 or hi > fa.get_reference_length(aln.reference_name):
            n_skip_ctx += 1
            continue
        ctx = fa.fetch(aln.reference_name, lo, hi).upper()

        if aln.is_reverse:
            read, ctxr = revcomp(seq), revcomp(ctx)
            clip5, clip3 = clipR, clipL
        else:
            read, ctxr = seq, ctx
            clip5, clip3 = clipL, clipR

        if low_complexity(read[clip5: L - clip3]):
            n_skip_lowcomp += 1
            continue

        # off = pad - clip5 is leakage-invariant: read pos p -> ctxr[off + p]
        off = pad - clip5
        rb = np.frombuffer(read.encode("ascii"), dtype=np.uint8)
        cb = np.frombuffer(ctxr.encode("ascii"), dtype=np.uint8)
        i = off + np.arange(L, dtype=np.int64)
        valid = (i >= 0) & (i < cb.size)
        matches = np.zeros(L, dtype=bool)
        matches[valid] = rb[valid] == cb[i[valid]]

        aln_s, aln_e = clip5, L - clip3
        # extend right through chance matches: a definition of the insert end
        # that depends on the genome, not on STAR's scoring
        e = aln_e
        while e < L and valid[e] and matches[e]:
            e += 1

        recs.append(_Read(seq=read, n=L, codes=_CODE[rb], matches=matches, valid=valid,
                          aln_s=aln_s, aln_e=aln_e, fp_end=e))
        anchor5 = (aln.reference_end + clip5) if aln.is_reverse \
            else (aln.reference_start - clip5)
        keys.append((aln.reference_name, aln.is_reverse, anchor5))
        loc_hist[keys[-1]] += 1

    recs, n_capped = cap_loci(recs, keys)
    return recs, loc_hist, n_seen, n_skip_cigar, n_skip_ctx, n_skip_lowcomp, n_capped


def cap_loci(recs: list[_Read], keys: list[tuple]) -> tuple[list[_Read], int]:
    """Keep at most a small, fixed number of reads per 5' locus.

    An abundant contaminant (an rRNA/tRNA/snoRNA fragment, an adapter dimer aligning to
    an adapter-like locus) is ONE molecule sequenced thousands of times. Left in, it
    turns the insert's own last bases into a "constant" block right in front of the
    adapter, and its genome context into the match profile. The read structure is a
    property of every molecule alike, so nothing is lost by letting each locus speak
    only a few times -- and no contaminant reference is needed to do it.
    """
    cap = max(LOCUS_CAP_MIN, int(LOCUS_CAP_FRAC * len(recs)))
    seen: Counter = Counter()
    kept = []
    for rec, key in zip(recs, keys):
        seen[key] += 1
        if seen[key] <= cap:
            kept.append(rec)
    return kept, len(recs) - len(kept)


def profile_window(recs: list[_Read], anchors: list[int | None], lo: int, hi: int):
    """Composition + genome-match rate at offsets [lo, hi) from a per-read anchor.

    Composition and match rate have different denominators: a read position
    always has a base, but it only has an implied genome base when the fetched
    context reached that far.
    """
    w = hi - lo
    comp = np.zeros((w, 4))
    match = np.zeros(w)
    n_comp = np.zeros(w)
    n_match = np.zeros(w)
    for rec, a in zip(recs, anchors):
        if a is None:
            continue
        base = a + lo                      # read position of window index 0
        j0 = max(0, -base)                 # ... clipped to the read
        j1 = min(w, rec.n - base)
        if j1 <= j0:
            continue
        sl = slice(base + j0, base + j1)
        n_comp[j0:j1] += 1
        codes = rec.codes[sl]
        for b in range(4):
            comp[j0:j1, b] += codes == b   # an N contributes to no base
        n_match[j0:j1] += rec.valid[sl]
        match[j0:j1] += rec.matches[sl]    # matches is False wherever valid is False
    return comp, match, n_comp, n_match


# --- de-novo adapter discovery ------------------------------------------------
def kmer_is_simple(km: str) -> bool:
    """Homopolymers, near-homopolymers and dinucleotide repeats cannot seed an
    adapter: they are tails (poly-A), dark cycles (poly-G) or simple repeats. Nor can
    a k-mer that is half tail (a run of >= 6): located with a partial match at the
    read end, it would anchor on every read that merely ends in the tail."""
    if max(km.count(b) for b in BASES) >= 0.75 * len(km) or len(set(km)) <= 2:
        return True
    return any(b * 6 in km for b in BASES)


def discover_seed(recs: list[_Read]):
    """The frequent complex k-mer lying NEAREST THE INSERT among the non-genomic 3'
    tails.

    Counted once per read. The 3' adapter is whatever is ligated nearest the insert,
    not whatever is most frequent: long reads run straight through a ligated linker
    into the sequencing-primer site, and then carry both in (nearly) every read. So
    all k-mers within `SEED_NEAR_FRAC` of the top count compete, and the one with the
    smallest mean distance from the insert end wins.

    -> (kmer | None, count, n_reads_with_room, top list of [kmer, count, mean_offset])
    """
    cnt: Counter = Counter()
    off: Counter = Counter()
    n_room = 0
    for rec in recs:
        if rec.n - rec.fp_end < SEED_K:
            continue
        n_room += 1
        read = rec.seq
        seen = set()
        for i in range(rec.fp_end, rec.n - SEED_K + 1):
            km = read[i: i + SEED_K]
            if km not in seen and "N" not in km:
                seen.add(km)
                cnt[km] += 1
                off[km] += i - rec.fp_end
    cands = [(k, c) for k, c in cnt.most_common(400) if not kmer_is_simple(k)]
    if not cands:
        return None, 0, n_room, []
    near = sorted((off[k] / c, k, c) for k, c in cands if c >= SEED_NEAR_FRAC * cands[0][1])
    _o, kmer, count = near[0]
    top = [[k, c, round(o, 2)] for o, k, c in near[:8]]
    return kmer, count, n_room, top


def homopolymer_tails(recs: list[_Read]) -> dict[str, int]:
    """Reads whose non-genomic tail opens with a homopolymer run, per base."""
    return {f"poly{b}": sum(1 for r in recs
                            if find_adapter(r.seq, b * HOMOPOLYMER_PROBE, r.fp_end - 2) >= 0)
            for b in BASES}


def anchor_on(recs: list[_Read], seq: str) -> list[int | None]:
    """Where `seq` starts in each read, or None.

    A short (>= 7 nt) partial match is only believed at or after the insert end -- a
    7-mer occurs by chance inside a footprint. A long match is specific wherever it
    starts, and it must be allowed to start *before* the computed insert end: when the
    genome happens to continue like the adapter (an abundant rRNA fragment ending in
    ...CTGTAGGCA, right in front of a linker that begins CTGTAGGCAC), the insert end is
    extended through those bases, and the true adapter start lies inside the "insert".
    """
    out: list[int | None] = []
    for r in recs:
        s = find_adapter(r.seq, seq, min_start=r.fp_end - SPECIFIC_OVERLAP,
                         min_overlap=SPECIFIC_OVERLAP)
        if s < 0:
            s = find_adapter(r.seq, seq, min_start=r.fp_end - 2)
        out.append(s if s >= 0 else None)
    return out


def gap_concentration(recs: list[_Read], anchors: list[int | None]) -> float:
    """How sharply the adapter sits at a fixed distance from the insert end.

    A real ligated adapter is separated from the insert by a construct of FIXED
    length (a UMI, a barcode, or nothing), so the gap has a sharp mode. A chance
    k-mer lands anywhere, so its gaps scatter.
    """
    gaps = Counter(a - r.fp_end for r, a in zip(recs, anchors) if a is not None)
    tot = sum(gaps.values())
    return max(gaps.values()) / tot if tot else 0.0


def assemble_scaffold(recs: list[_Read], anchors: list[int | None], seed: str) -> str:
    """Extend the seed to the full fixed 3' scaffold by per-position consensus.

    Left and right of the seed, a position belongs to the scaffold for as long as
    it is constant across the anchored reads. The left walk stops by itself at a
    UMI (random) or at the insert (genomic, so variable). A homopolymer run inside
    the left extension is a poly-A tail sitting between insert and adapter --
    constant in composition but variable in length -- so the scaffold starts after
    it, and sheds any leading bases of the tail's kind (whether the adapter itself
    begins with that base cannot be known, and does not matter for trimming).
    """
    comp, _m, n_comp, _n = profile_window(recs, anchors, -MAXW, MAXW)
    n_anch = sum(a is not None for a in anchors)
    frac = comp / np.maximum(comp.sum(axis=1, keepdims=True), 1)

    def const_base(j: int, min_n: float) -> str | None:
        if j < 0 or j >= len(frac) or n_comp[j] < min_n:
            return None
        if entropy(frac[j]) > SCAFFOLD_ENT:
            return None
        return BASES[int(frac[j].argmax())]

    left = ""
    for d in range(1, MAXW + 1):
        b = const_base(MAXW - d, max(30, 0.5 * n_anch))
        if b is None:
            break
        left = b + left
    right = ""
    for d in range(len(seed), MAXW):
        # far right of the anchor only the short-insert reads still have bases
        b = const_base(MAXW + d, max(30, 0.05 * n_anch))
        if b is None:
            break
        right += b

    # the last homopolymer run that BEGINS in the left extension or at the seed's first
    # base (the seed may open with a few tail bases) is the tail; runs deeper inside the
    # seed are the adapter's own sequence
    scaffold = left + seed + right
    tail_base = None
    run_end, i = -1, 0
    while i <= len(left) and i < len(scaffold):
        j = i
        while j < len(scaffold) and scaffold[j] == scaffold[i]:
            j += 1
        if j - i >= HOMOPOLYMER_MIN:
            run_end, tail_base = j, scaffold[i]
        i = j
    if run_end >= 0:
        scaffold = scaffold[run_end:]
    if tail_base:
        scaffold = scaffold.lstrip(tail_base)
    return scaffold


def rate(a: np.ndarray, b: np.ndarray) -> list:
    return np.where(b > 0, a / np.maximum(b, 1), np.nan).tolist()


def comp_frac(c: np.ndarray) -> list:
    s = c.sum(axis=1, keepdims=True)
    return np.where(s > 0, c / np.maximum(s, 1), np.nan).tolist()


def profile_bam(bam: str, fasta: str, *, label: str = "", max_reads: int = 200_000,
                min_anchor_frac: float = MIN_ANCHOR_FRAC) -> dict:
    """Positional match/composition profile of a locally-aligned BAM.

    Purely descriptive: every number here is a measurement, and the structure call
    (`fqdissect.infer`) is made from this dict alone, never from the BAM.
    """
    recs, loc_hist, n_seen, n_skip_cigar, n_skip_ctx, n_skip_lc, n_capped = load_reads(
        bam, fasta, max_reads)
    n = len(recs)
    if n == 0:
        raise ValueError(f"no usable alignments in {bam}")

    top5p_count = loc_hist.most_common(1)[0][1] if loc_hist else 0
    n_loaded = n + n_capped

    len_hist = Counter(r.n for r in recs)
    clip5_hist = Counter(r.aln_s for r in recs)
    clip3_hist = Counter(r.n - r.aln_e for r in recs)
    fp_len_hist = Counter(r.fp_end - r.aln_s for r in recs)
    tail_hist = Counter(r.n - r.fp_end for r in recs)   # non-genomic 3' length

    # ---- 5' side: anchor on read position 0 (exact, no leakage)
    c5, m5, _n5c, n5 = profile_window(recs, [0] * n, 0, MAXP5)

    # ---- 3' side, anchor A: the read's own 3' end (offset 0 = last base)
    ct, mt, _ntc, nt_ = profile_window(recs, [r.n - 1 for r in recs], -(MAXW - 1), 1)
    ct, mt, nt_ = ct[::-1], mt[::-1], nt_[::-1]   # index 0 = last base of read

    # ---- 3' side, anchor B: the adapter start, found de novo
    seed, seed_count, n_room, seed_top = discover_seed(recs)
    homo_hits = homopolymer_tails(recs)

    anchor_kind, scaffold = "none", ""
    anchors: list[int | None] = [None] * n
    seed_gap_conc = 0.0
    if seed and seed_count >= MIN_ANCHOR_READS:
        seed_anchors = anchor_on(recs, seed)
        # gate on the EXACT k-mer count, not on the located reads: locating allows a
        # partial match at the read end, which is generous to a seed that is not real
        hits = seed_count
        seed_gap_conc = gap_concentration(recs, seed_anchors)
        # A MINORITY of reads showing the adapter does not mean there is no adapter:
        # a long insert pushes it off the end of most reads. Accept it when there
        # are enough reads to measure a construct from and its distance to the
        # insert end is actually fixed.
        if hits >= SEED_MIN_FRAC * n or seed_gap_conc >= MIN_GAP_CONC:
            scaffold = assemble_scaffold(recs, seed_anchors, seed)
            anchors = anchor_on(recs, scaffold)
            # re-anchoring can only lose reads if the consensus went wrong
            if sum(a is not None for a in anchors) < 0.8 * sum(a is not None
                                                               for a in seed_anchors):
                LOG.warning("scaffold %s anchors fewer reads than its seed; keeping the seed",
                            scaffold)
                scaffold, anchors = seed, seed_anchors
            anchor_kind = "denovo"
    if anchor_kind == "none":
        # nothing but a homopolymer: report it so the caller can name the tail
        best = max(homo_hits, key=lambda k: homo_hits[k])
        if homo_hits[best] >= min_anchor_frac * n:
            anchor_kind, scaffold = f"homopolymer:{best}", best[-1] * HOMOPOLYMER_PROBE
            anchors = anchor_on(recs, scaffold)

    n_anchored = sum(a is not None for a in anchors)
    if n_anchored:
        ca, ma, nac, na = profile_window(recs, anchors, -MAXW, MAXW)
    else:
        ca = np.full((2 * MAXW, 4), np.nan)
        ma = np.full(2 * MAXW, np.nan)
        nac = na = np.zeros(2 * MAXW)

    # Distance from the insert end to the adapter start. A fixed-length construct
    # (UMI, barcode) gives a sharp mode; a poly(A) tail gives a broad distribution.
    gap_hist = Counter(a - r.fp_end for r, a in zip(recs, anchors) if a is not None)

    out = {
        "label": label or os.path.basename(bam),
        "n_alignments_seen": n_seen,
        "n_used": n,
        "n_skipped_indel_or_junction": n_skip_cigar,
        "n_skipped_contig_edge": n_skip_ctx,
        "n_skipped_low_complexity": n_skip_lc,
        "n_capped_abundant_locus": n_capped,
        # share of the loaded reads at the single busiest 5' locus, BEFORE capping
        "top5p_locus_frac": top5p_count / n_loaded,
        "n_distinct_5p_loci": len(loc_hist),
        "read_len_hist": dict(sorted(len_hist.items())),
        "clip5_hist": dict(sorted(clip5_hist.items())),
        "clip3_hist": dict(sorted(clip3_hist.items())),
        "footprint_len_hist": dict(sorted(fp_len_hist.items())),
        "tail_len_hist": dict(sorted(tail_hist.items())),
        "fpend_to_adapter_gap_hist": dict(sorted(gap_hist.items())),
        # 5'-anchored, index 0 = first base of the read
        "p5_comp": comp_frac(c5),
        "p5_match": rate(m5, n5),
        "p5_n": n5.tolist(),
        # read-3'-end-anchored, index 0 = last base of the read, growing inwards
        "t3_comp": comp_frac(ct),
        "t3_match": rate(mt, nt_),
        "t3_n": nt_.tolist(),
        # adapter-anchored, index MAXW = first base of the scaffold
        "anchor_kind": anchor_kind,
        "anchor_seq": scaffold,
        "anchor_offset": MAXW,
        "n_anchored": n_anchored,
        "frac_anchored": n_anchored / n,
        "adap_comp": comp_frac(ca),
        "adap_match": rate(ma, na),
        "adap_n": na.tolist(),
        "adap_n_comp": nac.tolist(),
        "seed_kmer": seed or "",
        "seed_count": seed_count,
        "seed_room": n_room,
        "seed_gap_concentration": seed_gap_conc,
        "seed_top": seed_top,
        "homopolymer_tail_hits": homo_hits,
        "bases": BASES,
    }
    LOG.info("%d usable alignments (%d more capped at abundant loci) | 3' anchor: %s %s "
             "(%.0f%% of reads) | skipped %d indel/junction, %d low-complexity", n, n_capped,
             anchor_kind, scaffold if anchor_kind == "denovo" else "",
             100 * out["frac_anchored"], n_skip_cigar, n_skip_lc)
    return out
