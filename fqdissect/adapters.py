"""Names for adapters -- used for LABELLING only, never for detection.

fqdissect finds the 3' adapter de novo (`fqdissect.profile`); nothing in this file
takes part in that. Once the fixed 3' scaffold has been assembled from the reads,
these cores are looked up *inside it* for two cosmetic-but-useful reasons:

* to put a familiar name on the adapter in the report, and
* to split a sample barcode off the front of the scaffold. A barcode is constant
  within one FASTQ, so from the reads alone `[barcode][adapter]` is one fixed
  block; only knowing where a standard adapter begins tells the two apart. Trimming
  is identical either way.

`--no-annotate` turns the lookup off; the call is then purely data-driven and the
whole scaffold is reported as the adapter.
"""
from __future__ import annotations

KNOWN_ADAPTERS = {
    "illumina_truseq": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "illumina_smallrna_ra3": "TGGAATTCTCGGGTGCCAAGG",
    "ingolia_linker": "CTGTAGGCACCATCAAT",
    "nextera": "CTGTCTCTTATACACATCT",
    "qiaseq_mirna": "AACTGTAGGCACCATCAAT",
}
CORE_LEN = 10
# constant blocks longer than this in front of a known adapter are a custom linker,
# not a sample barcode
MAX_BARCODE = 8

DISPLAY = {
    "illumina_truseq": "TruSeq",
    "illumina_smallrna_ra3": "smallRNA-RA3",
    "ingolia_linker": "Ingolia-linker",
    "nextera": "Nextera",
    "qiaseq_mirna": "QIAseq-miRNA",
    "denovo": "adapter",
    "none": "no-adapter",
    "unknown": "adapter?",
    "none_visible": "adapter beyond tail?",
}


def annotate(block: str) -> tuple[str, int] | None:
    """(name, offset) of the leftmost known adapter whose core occurs in `block`,
    or None -- also when it sits too deep in the block for the prefix to be a
    barcode (`CACTCGGGCACCAAGGAC` + TruSeq read-through is one 18 nt linker)."""
    best = None
    for name, seq in KNOWN_ADAPTERS.items():
        i = block.find(seq[:CORE_LEN])
        if i >= 0 and (best is None or i < best[1]):
            best = (name, i)
    if best and best[1] > MAX_BARCODE:
        return None
    return best


def display(name) -> str:
    return DISPLAY.get(str(name), str(name))
