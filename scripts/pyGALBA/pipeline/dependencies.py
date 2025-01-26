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
            missing_deps.append(f"Missing Python library '{lib}' (try 'pip install {lib}')")

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

    # We'll store all missing dependencies in this list of strings:
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

    # 3) Conditionals for RNA-seq / Iso-Seq:
    # If user provided FASTQ-based short reads (rnaseq_fq), we need HISAT2:
    if args.rnaseq_fq and len(args.rnaseq_fq) > 0:
        logger.info("User provided RNA-seq in FASTQ => Checking for 'hisat2' tools.")
        which_executable("hisat2", "HISAT2 (hisat2)", missing_deps)
        which_executable("hisat2-build", "HISAT2 (hisat2-build)", missing_deps)

    # If user provided FASTQ-based Iso-Seq => minimap2 is required:
    if args.isoseq_fq and len(args.isoseq_fq) > 0:
        logger.info("User provided Iso-Seq in FASTQ => Checking for 'minimap2'.")
        which_executable("minimap2", "minimap2", missing_deps)

    # If user provided any rnaseq_bam or isoseq_bam or rnaseq_fq or isoseq_fq, 
    # we might need samtools unless the user gave a StringTie GTF:
    need_samtools = False
    if (args.rnaseq_bam and len(args.rnaseq_bam) > 0) or \
       (args.isoseq_bam and len(args.isoseq_bam) > 0) or \
       (args.rnaseq_fq and len(args.rnaseq_fq) > 0) or \
       (args.isoseq_fq and len(args.isoseq_fq) > 0):
        # If user *didn't* provide a stringtie_gtf => we need samtools
        if not args.stringtie_gtf:
            need_samtools = True

    if need_samtools:
        logger.info("Transcript reads (RNA-seq or Iso-Seq) not in StringTie => Checking for 'samtools'.")
        which_executable("samtools", "samtools", missing_deps)

    # 4) For any RNASEQ or ISOSEQ that is not a pre-assembled stringtie -> we need stringtie:
    need_stringtie = False
    if not args.stringtie_gtf:
        # Means user must do assembly themselves => stringtie is needed
        if (args.rnaseq_fq and len(args.rnaseq_fq) > 0) or \
           (args.rnaseq_bam and len(args.rnaseq_bam) > 0) or \
           (args.isoseq_fq and len(args.isoseq_fq) > 0) or \
           (args.isoseq_bam and len(args.isoseq_bam) > 0):
            need_stringtie = True

    if need_stringtie:
        logger.info("No StringTie GTF => Need to run StringTie ourselves.")
        which_executable("stringtie", "StringTie", missing_deps)

    # 5) If any rnaseq or isoseq was provided in any format, we need TransDecoder + bedtools
    if (args.rnaseq_fq or args.rnaseq_bam or 
        args.isoseq_fq or args.isoseq_bam or 
        args.stringtie_gtf):
        logger.info("Transcript data was provided => Checking TransDecoder + bedtools.")
        which_executable("TransDecoder.LongOrfs", "TransDecoder.LongOrfs", missing_deps)
        which_executable("TransDecoder.Predict", "TransDecoder.Predict", missing_deps)
        # TransDecoder "util" scripts if needed:
        #   check "gtf_genome_to_cdna_fasta.pl" etc. if you call them by name
        which_executable("bedtools", "bedtools", missing_deps)

    # End of dependency checks: raise an error if anything is missing
    if missing_deps:
        # Format them into one multi-line or single-line error
        error_msg = (
            "The following dependencies are missing or not found:\n  - " 
            + "\n  - ".join(missing_deps)
            + "\nPlease install or configure them before running galba.py again."
        )
        # Could be a custom exception, or just sys.exit
        logger.error(error_msg)
        sys.exit(1)
    else:
        logger.info("All required tools appear to be available for the given pipeline configuration.")
