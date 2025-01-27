import logging
import os
import sys
import subprocess
import shutil

def run_miniprot(args):
    """
    Align proteins to the genome using Miniprot, produce an ALN,
    run boundary_scorer, then run miniprothint.py to get a final GTF.
    """

    logger = logging.getLogger("galba_pipeline.runminiprot")

    # 1) Prepare working directory
    workdir = args.workingdir
    os.makedirs(workdir, exist_ok=True)

    # 2) Build the Miniprot index
    index_file = os.path.join(workdir, "genome.mpi")
    build_cmd = [
        "miniprot",
        "-t", str(args.threads),
        "-d", index_file,
        args.genome
    ]
    logger.info("Building Miniprot index...")
    logger.debug(f"Command: {' '.join(build_cmd)}")
    try:
        subprocess.run(build_cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Miniprot index-building failed with exit code {e.returncode}.")
        sys.exit(1)

    # 3) Miniprot alignment
    aln_file = os.path.join(workdir, "miniprot.aln")
    align_cmd = [
        "miniprot",
        "-I",
        "-ut" + str(args.threads),
        "--aln",
        index_file,
        args.prot_seq
    ]
    logger.info("Running Miniprot alignment...")
    logger.debug(f"Command: {' '.join(align_cmd)}")

    try:
        with open(aln_file, "w") as out_aln:
            subprocess.run(align_cmd, stdout=out_aln, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Miniprot alignment failed with exit code {e.returncode}.")
        sys.exit(1)

    # 4) Boundary scoring with miniprot_boundary_scorer
    # Attempt to find the scorer in PATH or use a user-supplied path if needed
    boundary_scorer = getattr(args, "miniprot_boundary_scorer", None) or "miniprot_boundary_scorer"
    boundary_scorer = shutil.which(boundary_scorer) or boundary_scorer

    if not os.path.isfile(boundary_scorer):
        logger.error(
            f"Cannot locate 'miniprot_boundary_scorer' executable at {boundary_scorer}. "
            "Please make sure it is in your PATH or specify --miniprot_boundary_scorer."
        )
        sys.exit(1)

    # For the BLOSUM matrix:
    if getattr(args, "scoring_matrix", None):
        scoring_matrix = args.scoring_matrix
    else:
        # same directory as boundary_scorer
        scorer_dir = os.path.dirname(os.path.abspath(boundary_scorer))
        scoring_matrix = os.path.join(scorer_dir, "blosum62.csv")

    if not os.path.isfile(scoring_matrix):
        logger.error(f"Could not find scoring matrix file: {scoring_matrix}")
        sys.exit(1)

    scorer_output = os.path.join(workdir, "miniprot.gff")
    boundary_cmd = [
        boundary_scorer,
        "-o", scorer_output,
        "-s", scoring_matrix
    ]

    logger.info("Scoring Miniprot alignment boundaries...")
    logger.debug(f"Command: {' '.join(boundary_cmd)}")
    try:
        # Instead of shell=True + < file, we do:
        with open(aln_file, "rb") as aln_in:
            subprocess.run(boundary_cmd, check=True, stdin=aln_in)
    except subprocess.CalledProcessError as e:
        logger.error(f"Boundary scorer failed with exit code {e.returncode}.")
        sys.exit(1)

    # 5) miniprothint.py step
    miniprot_gtf = os.path.join(workdir, "miniprot.gtf")
    prothint_cmd = [
        "miniprothint.py",
        scorer_output,
        "--workdir", workdir,
        "--ignoreCoverage",
        "--topNperSeed", "10",
        "--minScoreFraction", "0.5"
    ]
    logger.info("Converting boundary-scored alignment to final GTF via miniprothint.py...")
    logger.debug(f"Command: {' '.join(prothint_cmd)}")

    try:
        subprocess.run(prothint_cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"miniprothint.py failed with exit code {e.returncode}.")
        sys.exit(1)

    # 6) Rename final GTF
    final_gtf_file = os.path.join(workdir, "miniprot_trainingGenes.gtf")
    produced_file = os.path.join(workdir, "miniprot.gtf")

    if os.path.isfile(produced_file):
        os.rename(produced_file, final_gtf_file)
        logger.info(f"Miniprot final GTF: {final_gtf_file}")
    else:
        logger.warning(
            "Expected 'miniprot.gtf' not found after miniprothint.py step."
        )

    logger.info("Miniprot alignment step complete.")
