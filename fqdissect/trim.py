"""Apply a read structure with cutadapt.

fqdissect does not trim reads itself. It translates the inferred structure into a
pipeline of cutadapt commands, records it as a shell script, and runs it. cutadapt
covers every structure fqdissect can call: it finds a (partial) adapter under an
explicit error rate, drops reads that lack it, trims a poly(A) tail exactly, and moves
a UMI from either end of the read into the header (`-u +-N --rename`).

What gets applied, in read order:

    5'  [linker][UMI]...            cut at fixed offsets; UMI blocks -> read header
        [barcode]                   found by anchored alignment (`-g ^BARCODE`, errors),
                                    so an indel in it moves the cut with it; reads
                                    without it are dropped
        [RT nt][== footprint ==]    kept
        [poly tail]                 trimmed (variable length)
        [nt][UMI][barcode]...       cut at fixed offsets from the scaffold; UMI -> header
        [barcode][adapter]          the scaffold: found by alignment, and removed with
                                    everything after it

The 3' blocks can only be located once the scaffold is gone, which is why a library
with a 3' UMI needs more than one pass (cutadapt applies `-u` *before* `-a`). The
passes are piped; nothing intermediate is written.

A 5' barcode is a fixed sequence, so it is located rather than counted: anything
fixed-length in front of it is cut in the same pass (cutadapt cuts `-u` before it
looks for adapters), anything behind it in the next one, counted from where the
barcode actually ended. Such a pass carries no 3' scaffold: cutadapt removes one
adapter per pass.

Reads that do not show the scaffold never reached the end of their molecule: their
3' end is the read's end, not the footprint's, and a 3' UMI beyond it was never
sequenced. They are dropped (`--discard-untrimmed`), which is standard practice for
Ribo-seq; `keep_untrimmed` is refused when a 3' UMI would come out short.
"""
from __future__ import annotations

import json
import os
import shlex
import subprocess

from . import structure
from .utils import LOG, ToolError, require_tools

ERROR_RATE = 0.12     # cutadapt's error model is a RATE: below ~9 nt of overlap the budget is 0
POLY_PROBE = 10       # poly-tail-only libraries: a run this long IS the 3' boundary
POLY_MIN = 6
# Errors allowed in a 5' barcode: 40% of its length (4 nt: 1, 5-7 nt: 2, 8 nt: 3).
# ERROR_RATE would allow none in a short barcode. On SRR24210493's 5-nt barcode (made
# antisense, docs/METHOD.md section 9), 1 error dropped 6k reads with two substitutions
# that a fixed cut handles; 2 errors keep them, for 140 more wrong cuts in 860k.
BARCODE5_ERROR_RATE = 0.4


def barcode_errors(seq: str) -> int:
    return max(1, int(BARCODE5_ERROR_RATE * len(seq)))


class UnsupportedStructure(ValueError):
    """The requested trimming cannot be done correctly for this read structure."""


# --- the plan ---------------------------------------------------------------------
def plan(call: dict) -> dict:
    """What has to happen to a read."""
    segs = structure.segments(call)
    # 5' steps in read order: fixed cuts (len, is_umi), merged, and anchored barcodes
    five: list = []
    run: list[tuple[int, bool]] = []
    for sg in segs:
        if sg.get("side") != 5 or sg["kind"] not in ("umi", "barcode", "linker"):
            continue
        if sg["kind"] == "barcode":
            five += _merge(run) + [("g", sg["seq"])]
            run = []
        else:
            run.append((sg["len"], sg["kind"] == "umi"))
    five += _merge(run)
    tail = call.get("polyA_tail", "none")
    poly = tail[-1] if tail not in ("none", "", None) else None
    # outermost first; `three` is walked from the read's 3' end inwards, once the
    # layer's scaffold has been removed
    layers = [{"scaffold": lay["scaffold"],
               "three": _merge([(s["len"], s["kind"] == "umi")
                                for s in reversed(lay["blocks"])])}
              for lay in structure.layers(call)]
    poly_is_anchor = False
    if layers[-1]["scaffold"] is None and poly and len(layers) == 1:
        # no adapter visible: the tail anchors itself. It sits between the footprint and
        # whatever follows, so cutting at it needs no knowledge of what follows.
        layers[-1]["scaffold"], poly_is_anchor = poly * POLY_PROBE, True
    return {"five": five, "layers": layers, "poly": poly, "poly_is_anchor": poly_is_anchor,
            "umi_len": sum(n for n, u in five if u is True)
                       + sum(n for lay in layers for n, u in lay["three"] if u)}


def _merge(ops: list[tuple[int, bool]]) -> list[tuple[int, bool]]:
    """Fuse neighbouring discards; UMI blocks stay separate from discards."""
    out: list[tuple[int, bool]] = []
    for n, umi in ops:
        if not n:
            continue
        if out and out[-1][1] == umi:
            out[-1] = (out[-1][0] + n, umi)
        else:
            out.append((n, umi))
    return out


def _five_passes(five: list) -> list[dict]:
    """Pack the 5' steps into passes: [one fixed cut]?[one anchored barcode]?, the cut
    first -- the order cutadapt applies them in."""
    passes: list[dict] = []
    for step in five:
        if step[0] == "g":
            if passes and "g" not in passes[-1]:
                passes[-1]["g"] = step[1]
            else:
                passes.append({"g": step[1]})
        else:
            passes.append({"five": step})
    return passes


def _three_passes(layers: list[dict]) -> list[dict]:
    """Pack the 3' work into cutadapt passes. Within one pass cutadapt cuts (`-u`) first
    and looks for the adapter (`-a`) second, so a pass is [one 3' cut]?[one scaffold]?."""
    events = []
    for lay in layers:
        if lay["scaffold"]:
            events.append(("a", lay["scaffold"]))
        events += [("u", op) for op in lay["three"]]
    passes: list[dict] = []
    for kind, val in events:
        if kind == "a" and passes and "a" not in passes[-1]:
            passes[-1]["a"] = val
        else:
            passes.append({kind: val})
    return passes


# --- command builders ---------------------------------------------------------------
def _cutadapt_passes(p: dict, *, in_fq: str, out_fq: str, min_len: int, threads: int,
                     min_overlap: int, discard_untrimmed: bool, log_base: str,
                     five: list, layers: list[dict]) -> list[dict]:
    """cutadapt passes: the 5' cuts, the 3' layers, and on the last one the poly tail
    and the minimum length."""
    # the two ends are independent; zip their passes, except that a pass with a 5'
    # barcode cannot also take a 3' scaffold (one adapter per read and pass)
    q5, q3 = _five_passes(five), _three_passes(layers)
    specs: list[dict] = []
    while q5 or q3 or not specs:
        spec = q5.pop(0) if q5 else {}
        if q3 and not ("g" in spec and "a" in q3[0]):
            spec.update(q3.pop(0))
        specs.append(spec)
    generic_poly = p["poly"] and not p["poly_is_anchor"] and p["poly"] != "A"
    if generic_poly and ("a" in specs[-1] or "g" in specs[-1]):
        specs.append({})                 # the generic poly trimmer is itself an `-a`
    have_umi = False

    stages = []
    for i, spec in enumerate(specs):
        cmd = ["cutadapt", "-j", str(threads)]
        umi_fields = []
        if spec.get("five"):
            n, is_umi = spec["five"]
            cmd += ["-u", str(n)]
            if is_umi:
                umi_fields.append("{cut_prefix}")
        if "u" in spec:
            n, is_umi = spec["u"]
            cmd += ["-u", str(-n)]
            if is_umi:
                umi_fields.append("{cut_suffix}")
        if "g" in spec:
            # reads without the barcode are dropped whatever --keep-untrimmed says: their
            # 5' cut would land at an unknown place
            cmd += ["-g", "^" + spec["g"], "-e", str(barcode_errors(spec["g"])),
                    "--discard-untrimmed"]
        if "a" in spec:
            mo = POLY_MIN if p["poly_is_anchor"] else min_overlap
            cmd += ["-a", spec["a"], "-O", str(mo), "-e", str(ERROR_RATE)]
            if discard_untrimmed:
                cmd.append("--discard-untrimmed")
        if umi_fields:
            cmd += ["--rename", "{id}" + ("" if have_umi else "_") + "".join(umi_fields)]
            have_umi = True
        last = i == len(specs) - 1
        if last:
            cmd += ["-m", str(min_len)]
            if p["poly"] and not p["poly_is_anchor"]:
                # cutadapt trims the poly tail AFTER -u and -a, which is the read order
                cmd += (["--poly-a"] if p["poly"] == "A" else
                        ["-a", p["poly"] * 100, "-O", str(POLY_MIN), "-e", "0.1"])
        log = f"{log_base}.cutadapt.{i + 1}"
        cmd += ["--json", log + ".json", "-o", out_fq if last else "-"]
        cmd.append("-" if i else in_fq)
        stages.append({"cmd": cmd, "log": log + ".log", "json": log + ".json"})
    return stages


def build_pipeline(call: dict, in_fq: str, out_fq: str, *, min_len: int = 20,
                   threads: int = 4, min_overlap: int = 7, discard_untrimmed: bool = True,
                   log_base: str | None = None) -> tuple[list[dict], dict]:
    """The piped cutadapt commands that apply `call` to `in_fq`. -> (stages, plan)"""
    if not call.get("functional"):
        raise ValueError(f"nothing to trim: the call is {call.get('status')!r} "
                         f"({call.get('reason', 'no structure')})")
    p = plan(call)
    if log_base is None:
        log_base = out_fq[:-3] if out_fq.endswith(".gz") else out_fq
        log_base = os.path.splitext(log_base)[0]
    if not discard_untrimmed and any(u for lay in p["layers"] for _, u in lay["three"]):
        raise UnsupportedStructure(
            "this library has a 3' UMI: a read without the adapter never sequenced it, so "
            "keeping untrimmed reads would write UMIs of differing length")
    return _cutadapt_passes(p, five=p["five"], layers=p["layers"], in_fq=in_fq,
                            out_fq=out_fq, min_len=min_len, threads=threads,
                            min_overlap=min_overlap, discard_untrimmed=discard_untrimmed,
                            log_base=log_base), p


def shell_script(stages: list[dict]) -> str:
    lines = ["#!/usr/bin/env bash", "# generated by fqdissect -- the exact pipeline that was run",
             "set -euo pipefail", ""]
    lines.append(" \\\n  | ".join(f"{shlex.join(s['cmd'])} 2> {shlex.quote(s['log'])}"
                                  for s in stages))
    return "\n".join(lines) + "\n"


# --- running ------------------------------------------------------------------------
def run_pipeline(stages: list[dict]) -> None:
    require_tools(*{s["cmd"][0] for s in stages})
    procs, logs = [], []
    prev = None
    try:
        for i, st in enumerate(stages):
            LOG.info("  %s%s", "| " if i else "$ ", shlex.join(st["cmd"]))
            log = open(st["log"], "w")
            logs.append(log)
            last = i == len(stages) - 1
            proc = subprocess.Popen(st["cmd"], stdin=prev,
                                    stdout=subprocess.DEVNULL if last else subprocess.PIPE,
                                    stderr=log)
            if prev is not None:
                prev.close()            # so an upstream tool gets SIGPIPE if we die
            prev = proc.stdout
            procs.append(proc)
        codes = [pr.wait() for pr in procs]
    finally:
        for log in logs:
            log.close()
        for pr in procs:
            if pr.poll() is None:
                pr.kill()
    for st, code in zip(stages, codes):
        if code:
            raise ToolError(f"{st['cmd'][0]} exited with status {code} (see {st['log']})")


def _counts(stage: dict) -> tuple[int | None, int | None]:
    try:
        with open(stage["json"]) as fh:
            rc = json.load(fh).get("read_counts", {})
    except (OSError, ValueError):
        return None, None
    return rc.get("input"), rc.get("output")


def trim_fastq(call: dict, in_fq: str, out_fq: str, **kw) -> dict:
    """Build, record and run the trimming pipeline. Returns its statistics."""
    stages, p = build_pipeline(call, in_fq, out_fq, **kw)
    os.makedirs(os.path.dirname(os.path.abspath(out_fq)) or ".", exist_ok=True)
    script = os.path.splitext(out_fq[:-3] if out_fq.endswith(".gz") else out_fq)[0] + ".sh"
    with open(script, "w") as fh:
        fh.write(shell_script(stages))
    LOG.info("trimming with cutadapt (%d pass%s)", len(stages), "" if len(stages) == 1 else "es")
    run_pipeline(stages)
    n_in = _counts(stages[0])[0]
    n_out = _counts(stages[-1])[1]
    stats = {"tool": "cutadapt", "script": script, "output": out_fq,
             "scaffolds": [lay["scaffold"] for lay in p["layers"]], "umi_len": p["umi_len"],
             "n_reads_in": n_in, "n_reads_out": n_out,
             "frac_kept": round(n_out / n_in, 4) if n_in and n_out is not None else None,
             "commands": [shlex.join(s["cmd"]) for s in stages]}
    if stats["frac_kept"] is not None:
        LOG.info("kept %s of %s reads (%.1f%%) -> %s", f"{n_out:,}", f"{n_in:,}",
                 100 * stats["frac_kept"], out_fq)
        if stats["frac_kept"] < 0.3:
            LOG.warning("most reads were dropped. If the molecule is about as long as the "
                        "read, the 3' scaffold runs off the end of most reads and they "
                        "cannot be cut at the footprint boundary; the survivors then skew "
                        "SHORT. Lowering --min-overlap recovers depth at the cost of "
                        "specificity.")
    return stats
