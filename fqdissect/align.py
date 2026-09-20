"""The permissive LOCAL alignment the read structure is read from.

The reads are aligned **untrimmed**: the adapter is an output of fqdissect, not an
input, so nothing may be removed before we look. STAR runs with
`--alignEndsType Local` and permissive length/score filters, so everything
non-genomic (5' UMI, RT additions, 3' UMI, barcode, adapter) is pushed into the
**soft clips** instead of preventing the alignment. These parameters are what make
the structure readable at all, so they are not exposed as options.
"""
from __future__ import annotations

import os
import shutil

from .utils import LOG, ToolError, require_tools, run

#   outFilterMultimapNmax 1       -- unique only; a multimapper has no single genomic
#                                    context to compare the read bases against
#   outFilterMatchNmin 20         -- 20 genomic bases is enough to place a footprint...
#   outFilterMatchNminOverLread 0 -- ...and the *fraction* filters must be off, or a
#   outFilterScoreMinOverLread 0     read that is half construct is thrown away for
#                                    being "too short" -- exactly the reads we need
#   outFilterMismatchNmax 3 / NoverLmax 0.12 -- tolerate real mismatches inside the
#                                    insert without letting construct bases in
#   seedSearchStartLmax 20        -- seed inside a short insert that is flanked by
#                                    non-genomic sequence
LOCAL_ARGS: list[str] = [
    "--outSAMtype", "BAM", "Unsorted",
    "--outSAMattributes", "NH", "HI", "AS", "nM", "MD",
    "--outSAMunmapped", "None",
    "--alignEndsType", "Local",
    "--outFilterMultimapNmax", "1",
    "--outFilterMatchNmin", "20",
    "--outFilterMatchNminOverLread", "0",
    "--outFilterScoreMinOverLread", "0",
    "--outFilterMismatchNmax", "3",
    "--outFilterMismatchNoverLmax", "0.12",
    "--seedSearchStartLmax", "20",
    "--alignSJoverhangMin", "8",
    "--alignSJDBoverhangMin", "2",
    "--outSJtype", "None",
]


def align_local(fastq: str, star_index: str, outdir: str, *, threads: int = 8,
                genome_load: str = "NoSharedMemory") -> str:
    """Align the (plain-text) sampled FASTQ; returns the unsorted BAM."""
    require_tools("STAR")
    if not os.path.isdir(star_index):
        raise ValueError(f"--star-index is not a directory: {star_index!r}")
    os.makedirs(outdir, exist_ok=True)
    # a killed STAR leaves _STARtmp behind and the next run refuses to start
    shutil.rmtree(os.path.join(outdir, "_STARtmp"), ignore_errors=True)
    cmd = ["STAR", "--runMode", "alignReads", "--runThreadN", str(int(threads)),
           "--genomeDir", star_index, "--genomeLoad", genome_load,
           "--readFilesIn", fastq,
           "--outFileNamePrefix", os.path.join(outdir, ""), *LOCAL_ARGS]
    LOG.info("STAR local alignment of the sample (loading the index can take minutes)")
    run(cmd, log_to=os.path.join(outdir, "star.log"))
    bam = os.path.join(outdir, "Aligned.out.bam")
    if not os.path.exists(bam) or os.path.getsize(bam) == 0:
        raise ToolError(f"STAR produced no alignments (see {outdir}/Log.out)")
    return bam
