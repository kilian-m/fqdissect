"""Synthetic libraries: a random "genome", inserts cut from it, a construct wrapped
around each insert, and a BAM that soft-clips the construct -- exactly what STAR's
local alignment hands to the profiler, without needing STAR or a real genome."""
from __future__ import annotations

import random

import pysam
import pytest

TRUSEQ = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
CUSTOM = "TCTCCTTGCATAATCACCAACC"          # in no table: must be found de novo
QIASEQ = "AACTGTAGGCACCATCAAT"


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_library(tmp_path, *, five="", three="", adapter=TRUSEQ, read_len=60, n=6000,
                 poly_a=False, downstream=None, seed=7, rt_frac=0.0):
    """`five`/`three` are templates: N = random base (UMI), ACGT = fixed (barcode).
    `downstream=(n_umi, second_const)` puts a UMI behind the adapter."""
    rng = random.Random(seed)
    genome = rand_seq(rng, 200_000)
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "\n".join(genome[i:i + 80] for i in range(0, len(genome), 80)) + "\n")
    pysam.faidx(str(fasta))

    def fill(t):
        return "".join(rng.choice("ACGT") if c == "N" else c for c in t)

    bam_path, fq_path = tmp_path / "lib.bam", tmp_path / "lib.fastq"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": len(genome)}]}
    truth = []
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam, open(fq_path, "w") as fq:
        for i in range(n):
            ins_len = rng.randint(26, 34)
            pos = rng.randint(1000, len(genome) - 1000)
            insert = genome[pos:pos + ins_len]
            f5, f3 = fill(five), fill(three)
            rt = rng.choice("AT") if rng.random() < rt_frac else ""
            tail = "A" * rng.randint(8, 16) if poly_a else ""
            down = (rand_seq(rng, downstream[0]) + downstream[1]) if downstream else ""
            full = f5 + rt + insert + tail + f3 + adapter + down + rand_seq(rng, read_len)
            read = full[:read_len]
            clip5 = len(f5) + len(rt)
            clip3 = read_len - clip5 - ins_len
            if clip3 < 0:
                continue
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = read
            a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            a.flag, a.reference_id, a.reference_start, a.mapping_quality = 0, 0, pos, 255
            a.cigar = [(op, ln) for op, ln in ((4, clip5), (0, ins_len), (4, clip3)) if ln]
            bam.write(a)
            fq.write(f"@r{i}\n{read}\n+\n{'I' * read_len}\n")
            umi = "".join(b for t, b in zip(five, f5) if t == "N") + \
                  "".join(b for t, b in zip(three, f3) if t == "N") + \
                  (down[:downstream[0]] if downstream else "")
            truth.append((f"r{i}", rt + insert, umi))
    return {"bam": str(bam_path), "fasta": str(fasta), "fastq": str(fq_path), "truth": truth}


@pytest.fixture
def library(tmp_path):
    def _make(**kw):
        return make_library(tmp_path, **kw)
    return _make
