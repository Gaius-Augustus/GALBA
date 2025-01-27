import shutil
import importlib
import sys
import logging

def check_python_libraries(missing_deps):
    """
    Make sure Python 3 + needed packages are installed: 
    BioPython, pandas, PyYAML, pygustus, etc.
    Instead of raising an exception immediately, we accumulate
    the missing dependencies in 'missing_deps' and return after checking all.
    """
    required_libraries = ["Bio", "pandas", "yaml", "pygustus"]
    for lib in required_libraries:
        spec = importlib.util.find_spec(lib)
        if spec is None:
            missing_deps.append(
                f"Missing Python library '{lib}' (try 'pip install {lib}')"
            )

def which_executable(exec_name, tool_label, missing_deps):
    """
    Returns the absolute path to 'exec_name' if found in PATH or None otherwise.
    If not found, add an entry to missing_deps with a message that 'tool_label' is missing.
    """
    path = shutil.which(exec_name)
    if path is None:
        missing_deps.append(
            f"Missing required tool '{tool_label}' (executable: '{exec_name}'). "
            "Please install or put into $PATH."
        )
    return path

def check_dependencies(args):
    """
    Inspect args to see which external tools are actually needed, 
    verify they are available, and list all missing items in one go.
    """

    logger = logging.getLogger("galba_pipeline.deps")
    missing_deps = []

    # 1) Always-required external tools:
    which_executable("diamond", "DIAMOND", missing_deps)
    which_executable("miniprot", "miniprot", missing_deps)
    which_executable("miniprot_boundary_scorer", "miniprot_boundary_scorer", missing_deps)
    which_executable("tsebra.py", "TSEBRA (tsebra.py)", missing_deps)

    # For AUGUSTUS tools:
    which_executable("augustus", "AUGUSTUS (augustus)", missing_deps)
    which_executable("etraining", "AUGUSTUS (etraining)", missing_deps)
    which_executable("optimize_augustus.pl", "optimize_augustus.pl", missing_deps)
    which_executable("new_species.pl", "new_species.pl", missing_deps)

    # 2) Required Python libraries:
    check_python_libraries(missing_deps)

    # 3) Conditionals for RNA-seq. 
    #    If user provided FASTQ-based short reads => we need HISAT2:
    needs_hisat = False
    if (args.rnaseq_single_fq and len(args.rnaseq_single_fq) > 0) \
       or (args.rnaseq_paired_fq and len(args.rnaseq_paired_fq) > 0):
        needs_hisat = True

    if needs_hisat:
        logger.info("User provided short reads => Checking for 'hisat2' tools.")
        which_executable("hisat2", "HISAT2 (hisat2)", missing_deps)
        which_executable("hisat2-build", "HISAT2 (hisat2-build)", missing_deps)

    # 4) If user provided Iso-Seq in FASTQ => we need minimap2:
    #    (If you use separate args for iso-seq, keep them here)
    if args.isoseq_fq and len(args.isoseq_fq) > 0:
        logger.info("User provided Iso-Seq in FASTQ => Checking for 'minimap2'.")
        which_executable("minimap2", "minimap2", missing_deps)

    # 5) If user provided any transcript reads (in FASTQ/BAM) but no pre-assembled GTF => 
    #    we need samtools (for BAM sorting/merging) and stringtie (assembly).
    #    For example:
    need_samtools = False
    need_stringtie = False

    # Check if there's any transcript data in any form:
    have_any_transcripts = (
        (args.rnaseq_single_fq and len(args.rnaseq_single_fq) > 0)
        or (args.rnaseq_paired_fq and len(args.rnaseq_paired_fq) > 0)
        or (args.rnaseq_bam and len(args.rnaseq_bam) > 0)
        or (args.isoseq_fq and len(args.isoseq_fq) > 0)
        or (args.isoseq_bam and len(args.isoseq_bam) > 0)
    )

    if have_any_transcripts and not args.stringtie_gtf:
        # => we are assembling ourselves => need samtools and stringtie
        need_samtools = True
        need_stringtie = True

    if need_samtools:
        logger.info("Transcript reads provided but no pre-assembled GTF => Checking samtools.")
        which_executable("samtools", "samtools", missing_deps)

    if need_stringtie:
        logger.info("No StringTie GTF => need to run StringTie ourselves.")
        which_executable("stringtie", "StringTie", missing_deps)

    # 6) If any transcripts (in any format), we need TransDecoder & bedtools
    if have_any_transcripts or args.stringtie_gtf:
        logger.info("Transcript data provided => Checking TransDecoder + bedtools.")
        which_executable("TransDecoder.LongOrfs", "TransDecoder.LongOrfs", missing_deps)
        which_executable("TransDecoder.Predict", "TransDecoder.Predict", missing_deps)
        which_executable("bedtools", "bedtools", missing_deps)
        # If you call TransDecoder utility scripts by name, check them too:
        # e.g. gtf_genome_to_cdna_fasta.pl, etc.

    # End of checks: if missing any dependencies, report them all at once:
    if missing_deps:
        error_msg = (
            "The following dependencies are missing or not found:\n  - " 
            + "\n  - ".join(missing_deps)
            + "\nPlease install or configure them before running galba.py again."
        )
        logger.error(error_msg)
        sys.exit(1)
    else:
        logger.info("All required tools appear to be available for the given pipeline configuration.")
