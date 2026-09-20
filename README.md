# fqdissect

Determine the **structure of the reads in a FASTQ file** — which part is the biological
insert, which parts are UMIs, barcodes, non-templated additions, poly(A) tail and
adapter — **de novo**, without being told the library protocol and without any table of
known adapter sequences.

```
$ fqdissect run SRR24210493.fastq.gz --genome GRCh38.fa --star-index STAR-index/ -o out --plot --trim
SRR24210493   ok   5'-[UMI,2nt]-[footprint,~30nt]-[UMI,5nt]-[barcode,TAGAC]-[TruSeq,AGATCGGAAGAG…]-3'
```

> **Scope (v0.1):** single-end **Ribo-seq** libraries. The method itself (below) is not
> Ribo-seq specific; the calling rules and their thresholds were tuned on Ribo-seq.

## How it works

1. **Sample** 200 k reads (uniform reservoir sample over the first 1 M).
2. **Align them untrimmed** with `STAR --alignEndsType Local` and permissive filters, so
   everything that is not genomic is pushed into the soft clips instead of preventing the
   alignment.
3. **Profile**: for every read position ask *does this base match the genome base the
   alignment implies for it?* Insert bases match ~95–99 %, construct bases ~25 % (chance).
   The 5′ side is anchored on read position 0; the 3′ side on the adapter start.
4. **Find the adapter de novo.** The alignment has already said where the insert ends,
   so k-mers are counted only in the *non-genomic 3′ tails* (k-mer counting in read tails
   is what fastp, DNApi, atropos `detect` and minion do blindly). Among the frequent
   k-mers the one **nearest the insert** seeds the adapter; its per-position consensus is
   extended left and right while positions stay constant. Homopolymers are tails, never
   adapters. Reads are capped per 5′ locus first, so an abundant rRNA/tRNA fragment cannot
   impersonate a constant block — no contaminant reference needed.
5. **Call** each position by three measurements — genome match (penetrance), base entropy
   and composition bias:

   | segment | signature | fate |
   |---|---|---|
   | footprint | genomic, variable | **keep** |
   | UMI | in every read, non-genomic, high entropy | trim → read header |
   | barcode | in every read, non-genomic, fixed base | trim |
   | RT untemplated nt | 1–2 nt, only in a *fraction* of reads, A/T-biased | **keep** |
   | template-switch / linker run | in every read, G/C-biased | trim |
   | poly(A) tail | homopolymer of variable length | trim |
   | adapter | constant block ending the molecule | trim |
   | UMI *behind* the adapter | random block between two constant blocks (QIAseq) | trim → read header |

   When the reads cannot support an answer the call is `undetermined` with a reason
   (exit status 3) — never a guess.

Details: [docs/METHOD.md](docs/METHOD.md).

## Install

```bash
conda env create -f environment.yml     # python, STAR, cutadapt, ...
conda activate fqdissect
pip install -e .
```

You need a genome FASTA and a STAR index of it (`STAR --runMode genomeGenerate`).

### External tools

| tool | needed for | why |
|---|---|---|
| **STAR** | the structure call (always) | Local alignment of *untrimmed* short reads with soft clipping — the whole method reads the structure out of the soft clips. |
| **cutadapt** | `--trim` only | It expresses every structure fqdissect can call: it finds a (partial) adapter under an explicit error rate, **drops reads that lack the adapter** (`--discard-untrimmed`), trims a poly(A) tail exactly (`--poly-a`), and moves a UMI from **either end** of the read into the header (`-u ±N --rename`). It is run, never imported. |

## Usage

```bash
fqdissect run reads.fastq.gz --genome genome.fa --star-index STAR-index/ -o out/ \
    [--plot] [--trim] [-p 16]

fqdissect plot out/reads.profile.json out/reads.structure.json      # redraw the figure
fqdissect trim reads.fastq.gz out/reads.structure.json -o trimmed.fastq.gz [--dry-run]
```

Useful options of `run`:

* `-p/--processes N` — handed to STAR, to pigz for decompression and to `cutadapt -j`.
  fqdissect's own profiling step takes a few seconds and stays serial; the run time is
  STAR loading its index plus the trimming of the full file.
* `-n/--sample`, `--scan` (0 = sample the whole file), `--bam` (reuse an existing local
  alignment), `--set NAME=VALUE` (override a calling threshold), `--no-annotate`,
  `--keep-tmp`.
* `--star-genome-load LoadAndKeep` shares the index between many runs; free it afterwards
  with `STAR --genomeLoad Remove --genomeDir …`.

### Outputs

| file | content |
|---|---|
| `<name>.structure.tsv` | the read structure as a table: one row per segment in read order, with length, fixed sequence and fate (see below) |
| `<name>.structure.json` | the full call: ordered `segments`, every measured quantity, thresholds used, flags, and the `adapter_scaffold` string to hand to a trimmer |
| `<name>.profile.json` | the positional statistics the call was made from |
| `<name>.structure.png` | `--plot`: the diagnostic figure |
| `<name>.trimmed.fastq.gz`, `<name>.trimmed.sh`, logs | `--trim` |
| stdout | `name <TAB> status <TAB> 5'-[…]-[…]-3'` |

```
# sample: SRR24210493
# status: ok
# structure: 5'-[UMI,2nt]-[footprint,~30nt]-[UMI,5nt]-[barcode,TAGAC]-[TruSeq,AGATCGGAAGAG…]-3'
order  segment    length  sequence              fate                       note
1      UMI        2                             trim, move to read header
2      footprint  20-40                         keep                       mode 30 nt
3      UMI        5                             trim, move to read header
4      barcode    5       TAGAC                 trim
5      adapter    40      AGATCGGAAGAGCACACG…   trim                       found de novo; matches TruSeq
```

(tab-separated in the file, sequences written out in full; the footprint length is the
observed 1st–99th percentile.)

### The figure

Every panel is a quantity a decision was made on, with the inferred segments and the
exact thresholds drawn on top — a visual audit of the call.

![example figure](docs/example_structure.png)

A: 5′ genome match and penetrance (the *step* is the 5′ UMI length) · B: 5′ composition
and entropy · C: 5′ soft-clip lengths · D: 3′ genome match, anchored on the de-novo
adapter · E: 3′ composition (the adapter consensus is read off here) · F: footprint
lengths.

### Trimming (`--trim`)

fqdissect does not trim reads itself: it writes a **cutadapt** pipeline, saves it as a
shell script (`<name>.trimmed.sh`, `--dry-run` prints it without running), and runs it:

```bash
cutadapt -j 16 -u 2 -a TAGACAGATCGGAAGAGC… -O 7 -e 0.12 --discard-untrimmed --rename '{id}_{cut_prefix}' -o - reads.fastq.gz \
  | cutadapt -j 16 -u -5 --rename '{id}{cut_suffix}' -m 20 -o out/reads.trimmed.fastq.gz -
```

More than one pass is needed whenever something sits at a fixed offset from the 3′
adapter (a 3′ UMI, a UMI behind the adapter): cutadapt cuts fixed lengths *before* it
looks for the adapter, so the 3′ cut has to wait for the pass after. The passes are
piped; nothing intermediate is written.

* UMIs end up as `@readname_UMI` (5′ blocks, then 3′ blocks) — the form `umi_tools dedup`
  and UMICollapse expect. The read comment is dropped.
* The untemplated RT base is kept with the footprint; reads shorter than `--min-len` are
  dropped.
* Reads that do not show the 3′ adapter never reached the end of their molecule — their
  3′ end is the read's end, and a 3′ UMI was never sequenced. They are dropped
  (`--discard-untrimmed`); `--keep-untrimmed` overrides this and is refused when the
  library has a 3′ UMI.
* A sample barcode is constant within one FASTQ, so from the data alone
  `[barcode][adapter]` is one fixed block. A small table of adapter *names*
  (`fqdissect/adapters.py`) is looked up **inside the de-novo block** to label the adapter
  and split the barcode off. It takes no part in detection; `--no-annotate` disables it.
  Trimming is identical either way.

## Tests

```bash
pytest          # synthetic libraries; needs no genome and no STAR
```

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).
