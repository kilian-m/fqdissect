"""Which way round the reads are: sense or antisense to the mRNA, or unstranded.

The reads alone cannot answer this. They show whether a library is stranded (each
locus is read from one strand), but a sense and an antisense library are both
strand-pure. The missing fact is which strand is transcribed, and it comes from an
annotation: each footprint inside a CDS is sense when it aligned to the CDS's
strand, antisense when it aligned to the other one.

Only run with `--gtf`. Without an annotation nothing is reported -- intron motifs
would do for eukaryotes, but prokaryotes have no introns, and a call that exists
for some genomes only is worse than none.

The unit is a distinct 5' locus, not a read: one abundant fragment over a CDS is one
molecule, and must not outvote the rest (cf. `profile.cap_loci`).
"""
from __future__ import annotations

import bisect
import gzip
from collections import defaultdict

import pysam

from .utils import LOG

UNSTRANDED_BAND = 0.10   # a sense share within 0.5 +/- this is an unstranded library


class Annotation:
    """Merged CDS intervals (or exons, for a GTF without CDS), per chromosome and strand."""

    def __init__(self, intervals: dict[tuple[str, str], list[tuple[int, int]]], feature: str,
                 path: str):
        self.feature = feature
        self.path = path
        self._iv: dict[tuple[str, str], tuple[list[int], list[int]]] = {}
        for key, ivs in intervals.items():
            starts, ends = [], []
            for s, e in sorted(ivs):
                if starts and s <= ends[-1]:
                    ends[-1] = max(ends[-1], e)
                else:
                    starts.append(s)
                    ends.append(e)
            self._iv[key] = (starts, ends)
        self.chroms = {c for c, _ in self._iv}
        self._alias: dict[str, str] = {}

    def match_references(self, references) -> int:
        """Map the genome's sequence names onto the GTF's (1 <-> chr1, MT <-> chrM).
        Returns how many genome sequences found a name in the GTF."""
        def candidates(name):
            yield name
            yield name[3:] if name.startswith("chr") else "chr" + name
            if name in ("MT", "chrM", "M", "chrMT"):
                yield from ("MT", "chrM", "M", "chrMT")
        self._alias = {}
        for ref in references:
            hit = next((c for c in candidates(ref) if c in self.chroms), None)
            if hit:
                self._alias[ref] = hit
        return len(self._alias)

    def strands_at(self, ref: str, blocks: list[tuple[int, int]]) -> set[str]:
        """The CDS strands overlapping any of the aligned blocks (0-based, half-open)."""
        chrom = self._alias.get(ref)
        found = set()
        if chrom is None:
            return found
        for strand in "+-":
            iv = self._iv.get((chrom, strand))
            if not iv:
                continue
            starts, ends = iv
            for s, e in blocks:
                i = bisect.bisect_left(starts, e) - 1
                if i >= 0 and ends[i] > s:
                    found.add(strand)
                    break
        return found


def load_gtf(path: str) -> Annotation:
    """CDS features of a GTF/GFF (plain or .gz); exons when it has no CDS at all."""
    feats: dict[str, dict] = {"CDS": defaultdict(list), "exon": defaultdict(list)}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 8 or f[2] not in feats or f[6] not in ("+", "-"):
                continue
            # GTF is 1-based inclusive -> 0-based half-open
            feats[f[2]][(f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
    feature = "CDS" if feats["CDS"] else "exon"
    if not feats[feature]:
        raise ValueError(f"--gtf {path}: no CDS or exon features")
    if feature == "exon":
        LOG.warning("--gtf has no CDS features; testing strandedness on exons")
    return Annotation(feats[feature], feature, path)


def call_strandedness(bam_path: str, ann: Annotation, *, min_loci: int,
                      strand_frac: float) -> dict:
    """Count footprints that lie in a CDS on one strand only, and call the library."""
    sense_loci, anti_loci = set(), set()
    n_reads = {"sense": 0, "antisense": 0}
    n_both = 0
    with pysam.AlignmentFile(bam_path) as bam:
        ann.match_references(bam.references)
        for aln in bam:
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if aln.has_tag("NH") and aln.get_tag("NH") > 1:
                continue
            hit = ann.strands_at(aln.reference_name, aln.get_blocks())
            if len(hit) != 1:
                n_both += len(hit) == 2
                continue
            read_strand = "-" if aln.is_reverse else "+"
            locus = (aln.reference_name, read_strand,
                     aln.reference_end if aln.is_reverse else aln.reference_start)
            if hit == {read_strand}:
                n_reads["sense"] += 1
                sense_loci.add(locus)
            else:
                n_reads["antisense"] += 1
                anti_loci.add(locus)

    n = len(sense_loci) + len(anti_loci)
    out = {"feature": ann.feature, "gtf": ann.path,
           "sense_loci": len(sense_loci), "antisense_loci": len(anti_loci),
           "sense_reads": n_reads["sense"], "antisense_reads": n_reads["antisense"],
           "reads_on_both_strands": n_both}
    if n < min_loci:
        out.update(call="undetermined", sense_fraction=(len(sense_loci) / n if n else None),
                   reason=f"only {n} footprint loci in a {ann.feature} on one strand "
                          f"(need {min_loci})")
        return out
    frac = len(sense_loci) / n
    out["sense_fraction"] = frac
    if frac >= strand_frac:
        out["call"] = "sense"
    elif frac <= 1 - strand_frac:
        out["call"] = "antisense"
    elif abs(frac - 0.5) <= UNSTRANDED_BAND:
        out["call"] = "unstranded"
    else:
        out["call"] = "undetermined"
        out["reason"] = (f"{frac:.0%} of {ann.feature} loci are sense: neither stranded "
                         f"nor unstranded")
    return out


def describe(st: dict) -> str:
    """One line for the table header and the log."""
    if st.get("call") == "undetermined":
        return f"undetermined ({st.get('reason', '')})"
    n = st["sense_loci"] + st["antisense_loci"]
    return f"{st['call']} ({st['sense_fraction']:.1%} of {n} {st['feature']} loci sense)"
