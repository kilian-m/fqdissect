"""The read structure as an ordered list of segments, and its renderings.

`segments(call)` is the single source of truth: the 5'->3' list of what a read is
made of. Everything a user or a downstream tool consumes is rendered from it:

* `to_tsv`            -- the segment table, one row per segment in read order;
* `structure_string`  -- the one-line human rendering;
* `fqdissect.trim`    -- the cutadapt commands.

Segment kinds and what happens to them:

    umi        random-templated   trimmed, moved to the read header (dedup)
    barcode    fixed-templated    trimmed
    linker     non-templated 5' run (template-switch G/C run, ...)   trimmed
    nontemplated   degenerate non-templated 3' bases                 trimmed
    rt         untemplated RT addition, in a fraction of reads       KEPT (part of the molecule)
    insert     the biological part (the ribosome footprint)          KEPT
    poly_tail  poly(A)-type tail, variable length                    trimmed
    adapter    the fixed 3' adapter, found de novo                   trimmed
"""
from __future__ import annotations

from . import adapters

_ABSENT = ("none", "unknown", "", None, "none_visible")


def _pct(hist: dict, q: float) -> int:
    items = sorted((int(k), v) for k, v in hist.items())
    tot = sum(v for _, v in items)
    acc = 0
    for k, v in items:
        acc += v
        if acc >= q * tot:
            return k
    return items[-1][0]


def segments(call: dict, profile: dict | None = None) -> list[dict]:
    """The read as ordered segments, 5'->3'. Empty for a call without a structure."""
    fn = call.get("functional")
    if not fn:
        return []
    segs: list[dict] = []

    # --- 5' construct: UMI / barcode blocks sit at known offsets; whatever else lies
    #     before the trim boundary is a non-templated run (template switch, linker)
    trim5 = int(fn.get("trim_5p", 0))
    pos = 0
    for blk in sorted(call.get("p5_layout") or [], key=lambda b: b["offset"]):
        if blk["offset"] > pos:
            segs.append({"kind": "linker", "side": 5, "len": blk["offset"] - pos})
        if blk["role"] == "barcode5":
            segs.append({"kind": "barcode", "side": 5, "len": blk["len"], "seq": blk["seq"]})
        else:
            segs.append({"kind": "umi", "side": 5, "len": blk["len"]})
        pos = blk["offset"] + blk["len"]
    if trim5 > pos:
        segs.append({"kind": "linker", "side": 5, "len": trim5 - pos})
    rt = int(fn.get("footprint_retains_rt_nt", 0))
    if rt:
        segs.append({"kind": "rt", "side": 5, "len": rt,
                     "penetrance": call.get("rt_penetrance")})

    # --- the insert
    hist = (profile or {}).get("footprint_len_hist") or {}
    ins = {"kind": "insert", "mode": call.get("footprint_len_mode")}
    if hist:
        ins["min_len"], ins["max_len"] = _pct(hist, 0.01), _pct(hist, 0.99)
    segs.append(ins)

    # --- 3' construct, footprint -> adapter
    tail = call.get("polyA_tail", "none")
    if tail not in _ABSENT:
        segs.append({"kind": "poly_tail", "side": 3, "base": tail[-1]})
    for blk in call.get("p3_layout") or []:
        kind = {"umi3": "umi", "barcode3": "barcode"}.get(blk["role"], "nontemplated")
        seg = {"kind": kind, "side": 3, "len": blk["len"]}
        if "seq" in blk:
            seg["seq"] = blk["seq"]
        segs.append(seg)
    if call.get("adapter3_name") not in _ABSENT and call.get("adapter3_seq") not in _ABSENT:
        segs.append({"kind": "adapter", "side": 3, "seq": call["adapter3_seq"],
                     "name": call["adapter3_name"]})
        n_down = int(call.get("umi3_downstream_len") or 0)
        if n_down:      # [adapter][UMI][second constant block] -- QIAseq-style
            seq2 = call["adapter3_downstream_seq"]
            hit = adapters.annotate(seq2) if call.get("adapter3_name") != "denovo" else None
            segs.append({"kind": "umi", "side": 3, "len": n_down, "downstream": True})
            segs.append({"kind": "adapter", "side": 3, "seq": seq2, "downstream": True,
                         "name": hit[0] if hit and hit[1] == 0 else "denovo"})
    return segs


# --- one-line rendering --------------------------------------------------------
def structure_string(call: dict, profile: dict | None = None) -> str:
    if call.get("status") != "ok":
        return f"{str(call.get('status', '?')).upper()}: {call.get('reason', '')}"
    parts = []
    for s in segments(call, profile):
        k = s["kind"]
        if k == "umi":
            parts.append(f"[UMI,{s['len']}nt]")
        elif k == "barcode":
            parts.append(f"[barcode,{s['seq']}]")
        elif k == "linker":
            parts.append(f"[linker,{s['len']}nt]")
        elif k == "nontemplated":
            parts.append(f"[nt,{s['len']}nt]")
        elif k == "rt":
            pen = s.get("penetrance")
            parts.append(f"[RT,{pen * 100:.0f}%]" if pen else f"[RT,{s['len']}nt]")
        elif k == "insert":
            parts.append(f"[footprint,~{s.get('mode')}nt]")
        elif k == "poly_tail":
            parts.append(f"[poly{s['base']}]")
        elif k == "adapter":
            seq = s["seq"] if len(s["seq"]) <= 13 else s["seq"][:12] + "…"
            parts.append(f"[{adapters.display(s['name'])},{seq}]")
    out = "5'-" + "-".join(parts) + "-3'"
    if call.get("deposit_state") not in ("raw", None):
        out += "   (3' adapter already trimmed)"
    return out


# --- what a trimmer needs ---------------------------------------------------------
def fixed_3p_blocks(call: dict) -> list[dict]:
    """The fixed-length 3' blocks still on the read once the scaffold is gone, in read
    order. A barcode directly in front of the adapter leaves with it (see
    `scaffold`), so it is not listed here."""
    blocks = [s for s in segments(call)
              if s.get("side") == 3 and s["kind"] in ("umi", "barcode", "nontemplated")
              and not s.get("downstream")]
    if blocks and blocks[-1]["kind"] == "barcode" and has_adapter(call):
        blocks = blocks[:-1]
    return blocks


def has_adapter(call: dict) -> bool:
    return (call.get("adapter3_name") not in _ABSENT
            and call.get("adapter3_seq") not in _ABSENT)


def scaffold(call: dict) -> str | None:
    """What an adapter trimmer should search for: the adapter, prefixed by the sample
    barcode when that sits directly in front of it. The barcode is just as constant,
    and it lies closer to the insert, so it is still visible in reads whose adapter
    has already run off the end."""
    if not has_adapter(call):
        return None
    lay = call.get("p3_layout") or []
    bc = lay[-1]["seq"] if lay and lay[-1]["role"] == "barcode3" else ""
    return bc + call["adapter3_seq"]


def layers(call: dict) -> list[dict]:
    """How the 3' end comes off, outermost first. Each layer is a scaffold to find by
    alignment (None: nothing to find, the blocks sit at the read end) and the
    fixed-length blocks that are exposed at the read's 3' end once it is gone.

    One layer, normally. Two when a UMI sits behind the adapter: first the constant
    block behind the UMI, which exposes the UMI; then the adapter itself.
    """
    inner = {"scaffold": scaffold(call), "blocks": fixed_3p_blocks(call)}
    down = [s for s in segments(call) if s.get("downstream")]
    if not down:
        return [inner]
    umi, second = down
    return [{"scaffold": second["seq"], "blocks": [umi]}, inner]


# --- the table -------------------------------------------------------------------
_FATE = {"umi": "trim, move to read header", "barcode": "trim", "linker": "trim",
         "nontemplated": "trim", "rt": "keep", "insert": "keep", "poly_tail": "trim",
         "adapter": "trim"}
_LABEL = {"umi": "UMI", "barcode": "barcode", "linker": "non-templated run",
          "nontemplated": "non-templated bases", "rt": "untemplated RT nt",
          "insert": "footprint", "poly_tail": "poly tail", "adapter": "adapter"}


def to_tsv(call: dict, profile: dict | None = None) -> str:
    """The read structure as a table, one segment per row in read order (5'->3').

    `length` is a number for fixed-length segments, `min-max` for the footprint (1st-99th
    percentile of what was observed) and for the RT addition (present in a fraction of
    the reads), and `variable` for a poly tail. `sequence` is given where it is fixed.
    """
    rows = []
    for i, s in enumerate(segments(call, profile), 1):
        k = s["kind"]
        if k == "insert":
            length = (f"{s['min_len']}-{s['max_len']}" if "min_len" in s else "variable")
            note = f"mode {s.get('mode')} nt"
        elif k == "rt":
            length = f"0-{s['len']}"
            pen = s.get("penetrance")
            note = f"in {pen * 100:.0f}% of reads" if pen else ""
        elif k == "poly_tail":
            length, note = "variable", ""
        elif k == "adapter":
            length = str(len(s["seq"]))
            note = "found de novo" + ("" if s["name"] == "denovo"
                                      else f"; matches {adapters.display(s['name'])}")
            if s.get("downstream"):
                note += "; constant block behind the UMI"
        else:
            length = str(s["len"])
            note = "behind the adapter" if s.get("downstream") else ""
        seq = s.get("seq") or (s["base"] * 3 + "..." if k == "poly_tail" else "")
        rows.append([str(i), _LABEL[k], length, seq, _FATE[k], note])
    head = [f"# sample: {call.get('sample', '')}",
            f"# status: {call.get('status', '')}"
            + (f" ({call['reason']})" if call.get("reason") else ""),
            f"# structure: {structure_string(call, profile)}",
            "\t".join(["order", "segment", "length", "sequence", "fate", "note"])]
    return "\n".join(head + ["\t".join(r) for r in rows]) + "\n"
