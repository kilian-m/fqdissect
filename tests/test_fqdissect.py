from __future__ import annotations

import gzip
import re
import shutil

import pytest

from conftest import CUSTOM, QIASEQ, TRUSEQ
from fqdissect import infer, profile, revcomp, strand, structure, trim
from fqdissect.cli import main
from fqdissect.utils import read_json

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


# --- strandedness (--gtf) ---------------------------------------------------------------
# every synthetic read aligns to the + strand of chr1, between 1000 and 199000
def _gtf(tmp_path, *rows, name="ann.gtf"):
    path = tmp_path / name
    path.write_text("".join(f"{c}\ttest\tCDS\t{s}\t{e}\t.\t{st}\t0\tgene_id \"g\";\n"
                            for c, s, e, st in rows))
    return str(path)


@pytest.mark.parametrize("rows, want", [
    ([("chr1", 1, 200_000, "+")], "sense"),
    ([("chr1", 1, 200_000, "-")], "antisense"),
    ([("chr1", 1, 100_000, "+"), ("chr1", 100_001, 200_000, "-")], "unstranded"),
    ([("1", 1, 200_000, "+")], "sense"),          # Ensembl names on a UCSC-named genome
])
def test_strandedness_from_the_annotated_cds_strand(library, tmp_path, rows, want):
    lib = library()
    ann = strand.load_gtf(_gtf(tmp_path, *rows))
    st = strand.call_strandedness(lib["bam"], ann, min_loci=200, strand_frac=0.8)
    assert st["call"] == want and st["feature"] == "CDS"
    assert st["sense_loci"] + st["antisense_loci"] > 5000


def test_strandedness_refuses_on_too_few_loci(library, tmp_path):
    lib = library()
    ann = strand.load_gtf(_gtf(tmp_path, ("chr1", 1, 3000, "+")))
    st = strand.call_strandedness(lib["bam"], ann, min_loci=200, strand_frac=0.8)
    assert st["call"] == "undetermined" and "loci" in st["reason"]


def test_strandedness_is_only_reported_with_a_gtf(library, tmp_path):
    lib = library(five="NN", three="NNNNNTAGAC")
    common = ["run", lib["fastq"], "--genome", lib["fasta"], "--bam", lib["bam"],
              "--set", "min_reads=1000", "-p", "2"]
    assert main(common + ["-o", str(tmp_path / "a"), "--name", "s"]) == 0
    assert "strandedness" not in read_json(tmp_path / "a" / "s.structure.json")
    gtf = _gtf(tmp_path, ("chr1", 1, 200_000, "-"))
    assert main(common + ["-o", str(tmp_path / "b"), "--name", "s", "--gtf", gtf]) == 0
    assert read_json(tmp_path / "b" / "s.structure.json")["strandedness"]["call"] == "antisense"
    assert "# strandedness: antisense" in (tmp_path / "b" / "s.structure.tsv").read_text()
    # a GTF for another genome fails before any work is done
    other = _gtf(tmp_path, ("scaffold_7", 1, 1000, "+"), name="other.gtf")
    assert main(common + ["-o", str(tmp_path / "c"), "--gtf", other]) == 2


# --- antisense libraries: the RNA's 5' end is at the read's 3' side -------------------
def test_antisense_rt_base_is_found_between_footprint_and_3p_umi(library):
    lib = library(five="NNNN", three="NNNN", adapter=CUSTOM, rt3_frac=0.9)
    prof = profile.profile_bam(lib["bam"], lib["fasta"], label="lib")
    call = infer.infer(prof, THR, antisense=True)
    assert call["status"] == "ok" and call["read_orientation"] == "antisense"
    assert call["rt3_len"] == 1 and call["umi3_len"] == 4 and call["umi5_len"] == 4
    assert call["functional"]["dedup_umi_len"] == 8
    assert re.search(r"-\[footprint,~\d+nt\]-\[RT,\d+%\]-\[UMI,4nt\]-",
                     structure.structure_string(call, prof))
    # the sense rules read the same base as a fifth UMI position, and trim it
    sense = infer.infer(prof, THR)
    assert sense["rt3_len"] == 0 and sense["umi3_len"] == 5


def test_antisense_rules_do_not_invent_an_rt_base(library):
    _, call = dissect(library(five="NN", three="NNNNNTAGAC"), antisense=True)
    assert call["status"] == "ok"
    assert call["rt3_len"] == 0 and call["umi3_len"] == 5 and call["umi5_len"] == 2


def test_antisense_poly_t_start_is_refused(library):
    _, call = dissect(library(poly_t5=True), antisense=True)
    assert call["status"] == "undetermined" and "antisense_5p_poly_tail" in call["reason"]


@pytest.mark.parametrize("keep", [False, True])
def test_antisense_trimming_writes_sense_footprints(library, tmp_path, keep):
    if not (shutil.which("cutadapt") and shutil.which("pigz")):
        pytest.skip("cutadapt or pigz not installed")
    lib = library(n=3000, five="NNNN", three="NNNN", adapter=CUSTOM, rt3_frac=0.9)
    _, call = dissect(lib, antisense=True)
    out = str(tmp_path / "trimmed.fastq.gz")
    stats = trim.trim_fastq(call, lib["fastq"], out, threads=2, keep_orientation=keep)
    assert stats["reverse_complemented"] is not keep and stats["n_reads_out"] > 2900
    got = _read_fastq(out)
    # the RT base stays on the footprint; reverse-complemented, it is back at its 5' end
    orient = (lambda x: x) if keep else profile.revcomp
    want = {f"{name}_{umi}": orient(insert) for name, insert, umi in lib["truth"]}
    assert sum(want.get(k) == v for k, v in got.items()) / len(got) > 0.97


def test_only_antisense_libraries_are_reverse_complemented(library):
    _, sense = dissect(library(five="NN", three="NNNNNTAGAC"))
    assert not any("fqdissect.revcomp" in c for c in _commands(sense))
    _, anti = dissect(library(five="NN", three="NNNNNTAGAC"), antisense=True)
    stages, _ = trim.build_pipeline(anti, "in.fq.gz", "out.fq.gz", threads=3)
    assert stages[-3]["cmd"][-3:-1] == ["-o", "-"]          # [..., -o, -, input]
    assert stages[-2]["cmd"][-2:] == ["-m", "fqdissect.revcomp"]
    assert stages[-1]["cmd"] == ["pigz", "-p", "3", "-1", "-c"]
    assert stages[-1]["stdout"] == "out.fq.gz"
    assert trim.shell_script(stages).rstrip().endswith("> out.fq.gz")


@pytest.mark.parametrize("chunk", [7, 64, 1 << 20])
def test_revcomp_handles_any_chunk_boundary(chunk):
    import io
    fq = b"".join(b"@r%d x\nACGTNacgtR%s\n+\nABCDEFGHIJ%s\n" % (i, b"A" * i, b"I" * i)
                  for i in range(50))
    out = io.BytesIO()
    assert revcomp.revcomp_fastq(io.BytesIO(fq), out, chunk=chunk) == 50
    recs = out.getvalue().split(b"\n")
    assert recs[1] == b"YacgtNACGT" and recs[3] == b"JIHGFEDCBA" and recs[0] == b"@r0 x"
    assert recs[5] == b"TYacgtNACGT"
    # a missing final newline is fine, a missing line is not
    assert revcomp.revcomp_fastq(io.BytesIO(fq.rstrip(b"\n")), io.BytesIO(), chunk) == 50
    with pytest.raises(ValueError):
        revcomp.revcomp_fastq(io.BytesIO(fq + b"@x\nACGT\n"), io.BytesIO(), chunk)


# --- a 5' barcode is located, not counted ---------------------------------------------
def test_a_5p_barcode_is_cut_by_anchored_alignment(library):
    _, call = dissect(library(five="NNGTCTANNNNN"))
    one, two = _commands(call)
    # the UMI in front is cut in the same pass, the one behind it in the next
    assert "-u 2 " in one and "-g ^GTCTA -e 2 --discard-untrimmed" in one
    assert "-a " not in one                       # one adapter per pass
    assert "-u 5 " in two and "-a AGATCGGAAGAG" in two


def test_barcode_errors_scale_with_its_length():
    assert [trim.barcode_errors("A" * n) for n in (3, 4, 5, 7, 8)] == [1, 1, 2, 2, 3]


def test_an_indel_in_a_5p_barcode_moves_the_cut(library, tmp_path):
    if not shutil.which("cutadapt"):
        pytest.skip("cutadapt not installed")
    lib = library(n=4000, five="GTCTANNNNN", five_del_frac=0.03)
    _, call = dissect(lib)
    assert call["status"] == "ok" and call["barcode5_seq"] == "GTCTA"
    out = str(tmp_path / "trimmed.fastq.gz")
    trim.trim_fastq(call, lib["fastq"], out, threads=2)
    got = _read_fastq(out)
    want = {f"{name}_{umi}": insert for name, insert, umi in lib["truth"]}
    # a fixed 10-nt cut gets every read with a deletion wrong (3%). Alignment fixes
    # most: a deleted LAST barcode base (GTCT+x) aligns as well as a substitution in
    # it, and is read as one -- the far commoner error -- so ~1/4 of deletions remain
    assert sum(want.get(k) == v for k, v in got.items()) / len(got) > 0.985
