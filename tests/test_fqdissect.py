from __future__ import annotations

import gzip
import shutil

import pytest

from conftest import CUSTOM, QIASEQ, TRUSEQ
from fqdissect import infer, profile, structure, trim
from fqdissect.cli import main

THR = infer.Thresholds(min_reads=1000)


def dissect(lib, **kw):
    prof = profile.profile_bam(lib["bam"], lib["fasta"], label="lib")
    return prof, infer.infer(prof, THR, **kw)


# --- de-novo adapter detection ------------------------------------------------------
def test_an_adapter_in_no_table_is_found_de_novo(library):
    prof, call = dissect(library(adapter=CUSTOM))
    assert prof["anchor_kind"] == "denovo"
    assert call["status"] == "ok"
    assert call["adapter3_name"] == "denovo"
    assert call["adapter3_seq"].startswith(CUSTOM)
    assert call["functional"]["trim_5p"] == 0 and call["umi3_len"] == 0


def test_the_known_adapter_table_only_labels(library):
    lib = library(three="NNNNNTAGAC")
    _, named = dissect(lib)
    _, blind = dissect(lib, annotate=False)
    # with the table: the barcode is split off a named adapter ...
    assert named["adapter3_name"] == "illumina_truseq"
    assert named["barcode3_seq"] == "TAGAC" and named["umi3_len"] == 5
    # ... without it the same bases are one de-novo scaffold; what gets trimmed is identical
    assert blind["adapter3_name"] == "denovo" and blind["umi3_len"] == 5
    assert blind["adapter3_seq"].startswith("TAGAC" + TRUSEQ[:12])
    assert structure.scaffold(named)[:20] == structure.scaffold(blind)[:20]


def test_umis_on_both_sides_and_a_barcode(library):
    prof, call = dissect(library(five="NN", three="NNNNNTAGAC"))
    assert call["umi5_len"] == 2 and call["umi3_len"] == 5
    assert call["functional"]["dedup_umi_len"] == 7
    assert structure.structure_string(call, prof).startswith(
        "5'-[UMI,2nt]-[footprint,~")


def test_the_adapter_nearest_the_insert_wins_over_read_through(library):
    # long reads run through the ligated linker into the sequencing-primer site, so
    # both are in every read; the 3' adapter is the one next to the insert
    _, call = dissect(library(adapter="CTGTAGGCACCATCAAT" + TRUSEQ, read_len=90))
    assert call["adapter3_name"] == "ingolia_linker"
    assert call["adapter3_seq"].startswith("CTGTAGGCACCATCAAT")


def test_a_umi_behind_the_adapter_is_found(library):
    prof, call = dissect(library(adapter=QIASEQ, downstream=(12, TRUSEQ), read_len=100))
    assert call["adapter3_seq"] == QIASEQ
    assert call["umi3_downstream_len"] == 12
    assert call["functional"]["dedup_umi_len"] == 12
    lays = structure.layers(call)
    assert len(lays) == 2 and lays[0]["scaffold"].startswith(TRUSEQ[:12])
    assert lays[1]["scaffold"] == QIASEQ


def test_a_poly_a_tail_is_a_tail_not_an_adapter(library):
    _, call = dissect(library(poly_a=True))
    assert call["polyA_tail"] == "polyA"
    assert call["adapter3_name"] == "illumina_truseq"
    assert not call["adapter3_seq"].startswith("AA")


def test_an_abundant_contaminant_cannot_fake_a_barcode(library, tmp_path):
    import pysam
    lib = library(adapter=CUSTOM)
    # one molecule, sequenced as often as everything else together
    src = pysam.AlignmentFile(lib["bam"])
    reads = list(src)
    out = str(tmp_path / "dominated.bam")
    with pysam.AlignmentFile(out, "wb", template=src) as fh:
        for r in reads:
            fh.write(r)
            fh.write(reads[0])
    prof = profile.profile_bam(out, lib["fasta"], label="x")
    call = infer.infer(prof, THR)
    assert prof["top5p_locus_frac"] > 0.45 and prof["n_capped_abundant_locus"] > 5000
    assert call["status"] == "ok" and call["barcode3_seq"] == "none"
    assert call["adapter3_seq"].startswith(CUSTOM)


def test_infer_refuses_rather_than_guesses():
    out = infer.infer({"label": "x", "n_used": 12})
    assert out["status"] == "undetermined" and "usable alignments" in out["reason"]
    with pytest.raises(ValueError):
        infer.infer({"n_used": 5})


def test_unknown_threshold_is_an_error():
    assert infer.Thresholds.from_overrides(["min_reads=10"]).min_reads == 10
    with pytest.raises(ValueError):
        infer.Thresholds.from_overrides(["no_such_knob=1"])


# --- the readable output ---------------------------------------------------------------
def test_the_table_lists_the_segments_in_read_order(library):
    prof, call = dissect(library(five="NN", three="NNNNNTAGAC"))
    lines = structure.to_tsv(call, prof).splitlines()
    assert lines[2].startswith("# structure: 5'-[UMI,2nt]-[footprint")
    rows = [ln.split("\t") for ln in lines if not ln.startswith("#")]
    assert rows[0] == ["order", "segment", "length", "sequence", "fate", "note"]
    assert [r[1] for r in rows[1:]] == ["UMI", "footprint", "UMI", "barcode", "adapter"]
    assert rows[1][2] == "2" and rows[1][4].endswith("read header")
    lo, hi = map(int, rows[2][2].split("-"))
    assert lo < hi and rows[2][4] == "keep"
    assert rows[4][3] == "TAGAC" and rows[5][3].startswith("AGATCGGAAGAG")


def test_process_count_reaches_the_tools(library):
    _, call = dissect(library(five="NNNN", adapter=CUSTOM))
    assert all("-j 12 " in c for c in _commands(call, threads=12))


# --- trimming with cutadapt -----------------------------------------------------------
def _commands(call, **kw):
    stages, _ = trim.build_pipeline(call, "in.fq.gz", "out.fq.gz", **kw)
    return [" ".join(s["cmd"]) for s in stages]


def test_cutadapt_cuts_the_3p_umi_only_after_the_scaffold_is_gone(library):
    _, call = dissect(library(five="NN", three="NNNNNTAGAC"))
    one, two = _commands(call)
    assert "-u 2 " in one and "-a TAGACAGATCGG" in one and "--discard-untrimmed" in one
    assert "{id}_{cut_prefix}" in one and "-u -5" not in one
    assert "-u -5" in two and "{id}{cut_suffix}" in two and "-m 20" in two


def test_a_poly_tail_anchors_the_cut_when_no_adapter_is_visible():
    call = {"status": "ok", "polyA_tail": "polyA", "adapter3_name": "none_visible",
            "adapter3_seq": "none_visible", "functional": {"trim_5p": 3}}
    (cmd,) = _commands(call)
    assert "-u 3 " in cmd and "-a AAAAAAAAAA -O 6" in cmd and "--poly-a" not in cmd


def test_keeping_untrimmed_reads_is_refused_when_a_3p_umi_would_come_out_short(library):
    _, call = dissect(library(three="NNNNN"))
    with pytest.raises(trim.UnsupportedStructure):
        _commands(call, discard_untrimmed=False)


def _read_fastq(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    return {lines[i][1:].split()[0]: lines[i + 1] for i in range(0, len(lines), 4)}


@pytest.mark.parametrize("lib_kw", [
    dict(five="NN", three="NNNNNTAGAC"),
    dict(adapter=QIASEQ, downstream=(12, TRUSEQ), read_len=100),
    dict(five="NNNN", adapter=CUSTOM, poly_a=True, read_len=80),
])
def test_trimmed_reads_are_the_inserts_and_the_header_carries_the_umi(library, tmp_path, lib_kw):
    if not shutil.which("cutadapt"):
        pytest.skip("cutadapt not installed")
    lib = library(n=3000, **lib_kw)
    _, call = dissect(lib)
    assert call["status"] == "ok"
    out = str(tmp_path / "trimmed.fastq.gz")
    stats = trim.trim_fastq(call, lib["fastq"], out, threads=2)
    got = _read_fastq(out)
    assert stats["n_reads_out"] == len(got) > 0.9 * len(lib["truth"])
    # an insert that itself ends in A is indistinguishable from the start of the tail
    end = (lambda x: x.rstrip("A")) if lib_kw.get("poly_a") else (lambda x: x)
    want = {f"{name}_{umi}": end(insert) for name, insert, umi in lib["truth"]}
    right = sum(want.get(name) == end(seq) for name, seq in got.items())
    assert right / len(got) > 0.97


# --- command line ---------------------------------------------------------------------
def test_run_from_a_bam_writes_every_output(library, tmp_path, capsys):
    lib = library(five="NN", three="NNNNNTAGAC")
    out = tmp_path / "out"
    rc = main(["run", lib["fastq"], "--genome", lib["fasta"], "--bam", lib["bam"],
               "-o", str(out), "--name", "s", "--plot", "--set", "min_reads=1000", "-p", "2",
               "--trim", "--dry-run"])
    assert rc == 0
    assert "5'-[UMI,2nt]-[footprint" in capsys.readouterr().out
    for ext in ("structure.json", "structure.tsv", "profile.json", "structure.png"):
        assert (out / f"s.{ext}").stat().st_size > 0
