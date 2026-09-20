"""Draw the read sample the structure is inferred from.

The read structure is a property of the library prep, identical in every read, so a
couple of hundred thousand reads determine it as well as the whole run does. The
sample is a uniform reservoir sample over the first `scan` reads (0 = the whole
file); archive FASTQs are in spot order, which is random with respect to content.
"""
from __future__ import annotations

import gzip
import random
import shutil
import subprocess

from .utils import LOG

# STAR is a short-read aligner whose input parser overflows (and segfaults) on
# multi-kilobase reads, so anything this long is not handed to it
MAX_READ_LEN = 500
MIN_READ_LEN = 18


def _open_stream(path: str, threads: int = 1):
    """(binary line stream, process-or-None). gzip is decompressed out of process,
    so stopping after `scan` reads does not pay for inflating the rest."""
    if path.endswith(".gz"):
        pigz = shutil.which("pigz")
        tool = pigz or shutil.which("gzip")
        if tool:
            cmd = [tool, "-dc"] + (["-p", str(max(1, threads))] if pigz else []) + [path]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL)
            return proc.stdout, proc
        return gzip.open(path, "rb"), None
    return open(path, "rb"), None


def sample_fastq(fastq: str, out_fastq: str, *, n: int = 200_000, scan: int = 1_000_000,
                 seed: int = 1, threads: int = 1) -> dict:
    """Reservoir-sample `n` reads from the first `scan` reads of `fastq` (0 = all)
    into the plain FASTQ `out_fastq`. Returns sampling statistics."""
    rng = random.Random(seed)
    limit = float("inf") if not scan else scan
    keep: list[tuple[bytes, bytes, bytes]] = []
    seen = n_bad = 0
    stream, proc = _open_stream(fastq, threads)
    try:
        while seen < limit:
            h, s, p, q = (stream.readline() for _ in range(4))
            if not h:
                break
            if not q or not h.startswith(b"@") or not p.startswith(b"+"):
                raise ValueError(f"{fastq}: malformed FASTQ record after {seen} reads")
            s, q = s.rstrip(), q.rstrip()
            if not (MIN_READ_LEN <= len(s) <= MAX_READ_LEN) or len(s) != len(q):
                n_bad += 1
                continue
            seen += 1
            rec = (h.split()[0], s, q)
            if len(keep) < n:
                keep.append(rec)
            else:
                j = rng.randrange(seen)
                if j < n:
                    keep[j] = rec
    finally:
        stream.close()
        if proc is not None:
            proc.kill()
            proc.wait()
    if not keep:
        raise ValueError(f"{fastq}: no usable reads ({n_bad} outside "
                         f"{MIN_READ_LEN}-{MAX_READ_LEN} nt)")
    rng.shuffle(keep)
    with open(out_fastq, "wb") as out:
        for h, s, q in keep:
            out.write(h + b"\n" + s + b"\n+\n" + q + b"\n")
    LOG.info("sampled %s of the first %s reads (%s skipped for length)",
             f"{len(keep):,}", f"{seen:,}", f"{n_bad:,}")
    return {"n_sampled": len(keep), "n_scanned": seen, "n_skipped_length": n_bad,
            "scan_limit": scan, "seed": seed}
