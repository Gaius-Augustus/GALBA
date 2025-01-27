# pipeline/runhisat.py

import logging
import os
import sys
import subprocess

def run_hisat(args):
    """
    Align RNA-Seq reads (FASTQ-based) to the genome using HISAT2,
    produce sorted BAM files, optionally merge them if single-end AND
    paired-end data are both present.
    """

    logger = logging.getLogger("galba_pipeline.runhisat")

    # If the user does NOT have short-read data, just log and return:
    # (Your pipeline might handle logic differently, e.g. skip or raise an error.)
    if not args.rnaseq_fq or len(args.rnaseq_fq) == 0:
        logger.info("No RNA-seq FASTQ data found => skipping HISAT2 alignment.")
        return

    # Possibly store index in a subfolder:
    index_prefix = os.path.join(workdir, "hisat_index", "genome")
    os.makedirs(os.path.join(workdir, "hisat_index"), exist_ok=True)

    # 2) Build the HISAT2 index
    #    If you have an environment variable or a path, you might do:
    #    hisat2_build = os.path.join(args.HISAT or "", "hisat2-build")
    #    but here we assume 'hisat2-build' is in PATH or found in check_dependencies.
    build_cmd = [
        "hisat2-build",
        "-p", str(args.threads),
        args.genome,
        index_prefix
    ]
    logger.info("Building HISAT2 index...")
    logger.debug(f"Command: {' '.join(build_cmd)}")

    try:
        subprocess.run(build_cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"HISAT2 index-building failed: {e}")
        sys.exit(1)

    # 3) Distinguish single-end vs paired-end read sets
    #    For simplicity, assume user can pass a list of read files (args.rnaseq_fq).
    #    If you'd like more advanced logic (split out pairs?), adapt as needed.
    #    We'll do a naive approach:
    single_files = []
    paired_files = []

    # Here we assume the user either gave single-end files or
    # gave them in pairs. We can parse that logic in arguments.py or here.
    # E.g. if user did: --rnaseq_fq file1.fq file2.fq ...
    # Possibly we do: single_files = [f for f in ... if ...], etc.

    # For a minimal example, let's pretend user must explicitly
    # provide the separate lists. If not, you can do some auto-detection.
    if getattr(args, "rnaseq_paired_sets", None):
        paired_files = args.rnaseq_paired_sets  # e.g. from your snippet
    if getattr(args, "rnaseq_single_sets", None):
        single_files = args.rnaseq_single_sets

    alignments_sam = []

    # 4) Align single-end short reads
    if single_files:
        logger.info("Mapping single-end RNA-seq reads with HISAT2...")
        # e.g. "hisat2 -x index_prefix -U file1,file2 --threads -S out.sam"
        single_out = os.path.join(workdir, "alignment_single_rnaseq.sam")
        single_files_str = ",".join(single_files)
        hisat_single_cmd = [
            "hisat2",
            "-x", index_prefix,
            "-U", single_files_str,
            "--dta",
            "-p", str(args.threads),
            "-S", single_out
        ]
        logger.debug(f"Command: {' '.join(hisat_single_cmd)}")
        try:
            subprocess.run(hisat_single_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"HISAT2 single-end mapping failed: {e}")
            sys.exit(1)
        alignments_sam.append(single_out)

    # 5) Align paired-end short reads
    if paired_files:
        # Typically, user passes them as [R1_1.fq, R1_2.fq, R2_1.fq, R2_2.fq, ...]
        # We'll assume they come in pairs. We'll do a simple approach:
        # Combine them all in one pass. If you'd prefer separate runs, adapt it:
        logger.info("Mapping paired-end RNA-seq reads with HISAT2...")
        paired_out = os.path.join(workdir, "alignment_paired_rnaseq.sam")

        # e.g. if user gave N files, we'd do all even as R1, all odd as R2
        # or in your snippet you do [0::2], [1::2].
        r1_list = paired_files[0::2]
        r2_list = paired_files[1::2]
        hisat_paired_cmd = [
            "hisat2",
            "-x", index_prefix,
            "-1", ",".join(r1_list),
            "-2", ",".join(r2_list),
            "--dta",
            "-p", str(args.threads),
            "-S", paired_out
        ]
        logger.debug(f"Command: {' '.join(hisat_paired_cmd)}")
        try:
            subprocess.run(hisat_paired_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"HISAT2 paired-end mapping failed: {e}")
            sys.exit(1)
        alignments_sam.append(paired_out)

    # 6) Convert all SAM => sorted BAM
    # We'll store them in the working dir as well:
    final_bams = []
    for samfile in alignments_sam:
        bamfile = os.path.splitext(samfile)[0] + ".bam"
        logger.info(f"Sorting SAM => BAM: {samfile} => {bamfile}")
        sort_cmd = [
            "samtools", "sort",
            "-@", str(args.threads),
            "-o", bamfile,
            samfile
        ]
        logger.debug(f"Command: {' '.join(sort_cmd)}")
        try:
            subprocess.run(sort_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"samtools sort failed: {e}")
            sys.exit(1)
        # delete the sam file
        os.remove(samfile)
        final_bams.append(bamfile)

    # Optionally, if we have multiple final_bams, merge them:
    if len(final_bams) > 1:
        logger.info("Merging multiple BAM files into single alignment_rnaseq.bam")
        merged_bam = os.path.join(workdir, "alignment_rnaseq.bam")
        merge_cmd = ["samtools", "merge", "-f", merged_bam] + final_bams
        logger.debug(f"Command: {' '.join(merge_cmd)}")
        try:
            subprocess.run(merge_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"samtools merge failed: {e}")
            sys.exit(1)

        # The final is alignment_rnaseq.bam
        # you can remove the per-file .bam if you like:
        for f in final_bams:
            os.remove(f)
    else:
        # If we only have one final bam, rename it to alignment_rnaseq.bam
        if final_bams:
            single_bam_path = final_bams[0]
            merged_bam = os.path.join(workdir, "alignment_rnaseq.bam")
            os.rename(single_bam_path, merged_bam)
        else:
            # means no reads => skip
            merged_bam = None

    logger.info("HISAT2 alignment step completed.")
