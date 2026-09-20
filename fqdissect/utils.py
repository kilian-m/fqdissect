"""Small shared helpers: logging, subprocesses, JSON, FASTQ streams."""
from __future__ import annotations

import gzip
import json
import logging
import os
import shlex
import shutil
import subprocess

import numpy as np

LOG = logging.getLogger("fqdissect")


def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(level=getattr(logging, level.upper(), logging.INFO),
                        format="%(asctime)s %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")


class ToolError(RuntimeError):
    """An external tool is missing or exited non-zero."""


def require_tools(*tools: str) -> None:
    missing = [t for t in tools if not shutil.which(t)]
    if missing:
        raise ToolError(f"required tool(s) not on PATH: {', '.join(missing)}")


def run(cmd: list[str], *, log_to: str | None = None) -> None:
    """Run one command; stdout+stderr go to `log_to` (appended) when given."""
    LOG.debug("$ %s", shlex.join(cmd))
    if log_to:
        with open(log_to, "a") as fh:
            fh.write(f"$ {shlex.join(cmd)}\n")
            fh.flush()
            proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
        tail = ""
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        tail = (proc.stdout or "")[-2000:]
    if proc.returncode:
        where = f" (see {log_to})" if log_to else f"\n{tail}"
        raise ToolError(f"{cmd[0]} exited with status {proc.returncode}{where}")


def open_fastq(path: str, mode: str = "rt"):
    """Open a FASTQ, transparently gzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, mode, compresslevel=4) if "w" in mode else gzip.open(path, mode)
    return open(path, mode)


def _jsonable(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o).__name__}")


def write_json(path: str, obj) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=_jsonable)
        fh.write("\n")
    os.replace(tmp, path)
    return path


def read_json(path: str):
    with open(path) as fh:
        return json.load(fh)


def sample_name(fastq: str) -> str:
    base = os.path.basename(fastq)
    for ext in (".gz", ".fastq", ".fq"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base
