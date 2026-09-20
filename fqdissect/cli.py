"""fqdissect command line.

    fqdissect run   reads.fastq.gz --star-index IDX --genome GENOME.fa -o out/ [--plot] [--trim]
    fqdissect plot  out/sample.profile.json out/sample.structure.json
    fqdissect trim  reads.fastq.gz out/sample.structure.json -o trimmed.fastq.gz
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile

from . import __version__
from .utils import LOG, ToolError, read_json, sample_name, setup_logging, write_json

LIBRARY_TYPES = ("riboseq",)


def _add_trim_options(ap: argparse.ArgumentParser) -> None:
    g = ap.add_argument_group("trimming (cutadapt)")
    g.add_argument("--min-len", type=int, default=20,
                   help="drop reads shorter than this after trimming [20]")
    g.add_argument("--min-overlap", type=int, default=7,
                   help="nt of the 3' scaffold that must be visible to cut on it [7]")
    g.add_argument("--keep-untrimmed", action="store_true",
                   help="keep reads that do not show the 3' adapter (default: drop them; "
                        "their 3' end is the read's end, not the footprint's)")
    g.add_argument("--dry-run", action="store_true",
                   help="write the trimming script but do not run it")
    g.add_argument("--force", action="store_true",
                   help="trim even when the structure call is 'undetermined' (best effort)")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="fqdissect",
        description="Determine the structure of the reads in a FASTQ file (biological "
                    "insert, UMIs, barcodes, adapter) de novo.")
    ap.add_argument("--version", action="version", version=f"fqdissect {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="infer the read structure of a FASTQ",
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    run.add_argument("fastq", help="single-end FASTQ (plain or .gz)")
    run.add_argument("-o", "--outdir", default="fqdissect_out")
    run.add_argument("--name", help="sample name [FASTQ basename]")
    run.add_argument("--library-type", choices=LIBRARY_TYPES, default="riboseq",
                     help="what the biological insert is")
    ref = run.add_argument_group("reference")
    ref.add_argument("--genome", required=True, help="genome FASTA (indexed, or indexable)")
    ref.add_argument("--star-index", help="STAR index of that genome (required unless --bam)")
    ref.add_argument("--bam", help="skip sampling+alignment: a BAM of the UNTRIMMED reads "
                                   "aligned with fqdissect's local STAR settings")
    smp = run.add_argument_group("sampling")
    smp.add_argument("-n", "--sample", type=int, default=200_000, help="reads to sample")
    smp.add_argument("--scan", type=int, default=1_000_000,
                     help="sample uniformly from the first SCAN reads (0 = whole file)")
    smp.add_argument("--seed", type=int, default=1)
    inf = run.add_argument_group("calling")
    inf.add_argument("--no-annotate", action="store_true",
                     help="do not look the de-novo adapter up in the table of known names "
                          "(the lookup only labels it and splits a sample barcode off its "
                          "front; detection never uses it)")
    inf.add_argument("--set", action="append", metavar="NAME=VALUE", default=[],
                     help="override a calling threshold (repeatable)")
    out = run.add_argument_group("outputs")
    out.add_argument("--plot", action="store_true", help="draw the diagnostic figure")
    out.add_argument("--plot-format", choices=("png", "pdf", "svg"), default="png")
    out.add_argument("--trim", action="store_true",
                     help="apply the structure to the full FASTQ with cutadapt: trim, and "
                          "move UMIs into the read header")
    out.add_argument("--keep-tmp", action="store_true", help="keep the sample and its BAM")
    run.add_argument("-p", "--processes", "-t", "--threads", dest="threads", type=int,
                     default=8, metavar="N",
                     help="processes/threads for STAR, gzip decompression (pigz) and "
                          "cutadapt (-j)")
    run.add_argument("--star-genome-load", default="NoSharedMemory",
                     help="STAR --genomeLoad (LoadAndKeep shares the index across runs)")
    run.add_argument("--log-level", default="INFO")
    _add_trim_options(run)

    pl = sub.add_parser("plot", help="(re)draw the diagnostic figure")
    pl.add_argument("profile", help="*.profile.json")
    pl.add_argument("structure", help="*.structure.json")
    pl.add_argument("-o", "--output", help="figure path [next to the profile, .png]")
    pl.add_argument("--log-level", default="INFO")

    tr = sub.add_parser("trim", help="apply an inferred structure to a FASTQ")
    tr.add_argument("fastq")
    tr.add_argument("structure", help="*.structure.json from `fqdissect run`")
    tr.add_argument("-o", "--output", required=True, help="trimmed FASTQ (.gz ok)")
    tr.add_argument("-p", "--processes", "-t", "--threads", dest="threads", type=int,
                    default=8, metavar="N",
                    help="processes/threads for cutadapt (-j)")
    tr.add_argument("--log-level", default="INFO")
    _add_trim_options(tr)
    return ap


def _trim(args, call: dict, fastq: str, out_fq: str) -> dict | None:
    from . import trim

    if call.get("status") != "ok" and not args.force:
        LOG.warning("not trimming: the structure is %s (--force applies the best-effort "
                    "plan)", call.get("status"))
        return None
    kw = dict(min_len=args.min_len, threads=args.threads, min_overlap=args.min_overlap,
              discard_untrimmed=not args.keep_untrimmed)
    if args.dry_run:
        stages, _ = trim.build_pipeline(call, fastq, out_fq, **kw)
        sys.stdout.write(trim.shell_script(stages))
        return None
    return trim.trim_fastq(call, fastq, out_fq, **kw)


def cmd_run(args) -> int:
    from . import align, infer, profile, sample, structure

    if not args.bam and not args.star_index:
        raise ValueError("give --star-index (or a pre-aligned --bam)")
    if not os.path.exists(args.fastq):
        raise ValueError(f"no such FASTQ: {args.fastq}")
    thr = infer.Thresholds.from_overrides(args.set)
    name = args.name or sample_name(args.fastq)
    os.makedirs(args.outdir, exist_ok=True)
    base = os.path.join(args.outdir, name)

    tmp = tempfile.mkdtemp(prefix=f"{name}.tmp.", dir=args.outdir)
    try:
        sampling = None
        bam = args.bam
        if not bam:
            fq = os.path.join(tmp, "sample.fastq")
            sampling = sample.sample_fastq(args.fastq, fq, n=args.sample, scan=args.scan,
                                           seed=args.seed, threads=args.threads)
            bam = align.align_local(fq, args.star_index, os.path.join(tmp, "star"),
                                    threads=args.threads, genome_load=args.star_genome_load)
        prof = profile.profile_bam(bam, args.genome, label=name,
                                   min_anchor_frac=thr.min_anchor_frac)
    finally:
        if not args.keep_tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    prof["sampling"] = sampling

    call = infer.infer(prof, thr, annotate=not args.no_annotate)
    call["fqdissect_version"] = __version__
    call["library_type"] = args.library_type
    call["fastq"] = os.path.abspath(args.fastq)
    call["structure"] = structure.structure_string(call, prof)
    call["segments"] = structure.segments(call, prof)
    if call["segments"]:
        call["adapter_scaffold"] = structure.scaffold(call)
    write_json(base + ".profile.json", prof)

    with open(base + ".structure.tsv", "w") as fh:
        fh.write(structure.to_tsv(call, prof))
    if args.plot:
        from . import plot
        fig = plot.plot_structure(prof, call, f"{base}.structure.{args.plot_format}")
        LOG.info("wrote %s", fig)
    if args.trim or args.dry_run:
        stats = _trim(args, call, args.fastq, base + ".trimmed.fastq.gz")
        if stats:
            call["trimming"] = stats
    write_json(base + ".structure.json", call)

    LOG.info("wrote %s", base + ".structure.json")
    for flag in call.get("flags") or []:
        LOG.info("flag: %s", flag)
    print(f"{name}\t{call['status']}\t{call['structure']}")
    return 0 if call["status"] == "ok" else 3


def cmd_plot(args) -> int:
    from . import plot

    out = args.output or args.profile.replace(".profile.json", "") + ".structure.png"
    print(plot.plot_structure(read_json(args.profile), read_json(args.structure), out))
    return 0


def cmd_trim(args) -> int:
    call = read_json(args.structure)
    stats = _trim(args, call, args.fastq, args.output)
    if stats:
        base = args.output[:-3] if args.output.endswith(".gz") else args.output
        write_json(os.path.splitext(base)[0] + ".trim.json", stats)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.log_level)
    try:
        return {"run": cmd_run, "plot": cmd_plot, "trim": cmd_trim}[args.command](args)
    except (ValueError, ToolError, OSError) as exc:
        LOG.error("%s", exc)
        return 2


if __name__ == "__main__":
    sys.exit(main())
