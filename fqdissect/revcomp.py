"""Reverse-complement a FASTQ stream: stdin -> stdout.

The last stage of `fqdissect trim` for an antisense library, so that the trimmed
footprints come out sense to the mRNA, the way every Ribo-seq tool expects them.
Read names, and the UMIs in them, are left as they are.

Pure Python on purpose: it runs concurrently with cutadapt and keeps up with it, so a
faster tool (seqtk, seqkit) would not shorten the trimming.

    python -m fqdissect.revcomp < trimmed.fastq > sense.fastq
"""
from __future__ import annotations

import sys

COMPLEMENT = bytes.maketrans(b"ACGTRYKMBDHVNacgtrykmbdhvn", b"TGCAYRMKVHDBNtgcayrmkvhdbn")
CHUNK = 1 << 22


def revcomp_fastq(inp, out, chunk: int = CHUNK) -> int:
    """Reverse-complement every record of a 4-line FASTQ. Returns the record count."""
    rest = b""
    n_records = 0
    eof = False
    while not eof:
        buf = inp.read(chunk)
        if not buf:
            eof = True
            if not rest.strip():
                break
            buf = b"\n"                          # a last line without its newline
        lines = (rest + buf).split(b"\n")
        n = (len(lines) - 1) // 4 * 4          # complete records; the tail waits
        rest = b"\n".join(lines[n:])
        if not n:
            continue
        lines = lines[:n]
        if not lines[0].startswith(b"@") or not lines[2].startswith(b"+"):
            raise ValueError("not a 4-line FASTQ")
        lines[1::4] = [s.translate(COMPLEMENT)[::-1] for s in lines[1::4]]
        lines[3::4] = [q[::-1] for q in lines[3::4]]
        out.write(b"\n".join(lines) + b"\n")
        n_records += n // 4
    if rest.strip():
        raise ValueError("truncated FASTQ: the last record is incomplete")
    return n_records


def main() -> int:
    try:
        revcomp_fastq(sys.stdin.buffer, sys.stdout.buffer)
    except ValueError as exc:
        print(f"fqdissect.revcomp: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
