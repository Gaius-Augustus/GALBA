import shutil
import importlib
import sys
import logging

def check_python_libraries():
    """
    Make sure Python 3 + needed packages are installed: 
    BioPython, pandas, PyYAML, pygustus, etc.
    """
    required_libraries = ["Bio", "pandas", "yaml", "pygustus"]
    for lib in required_libraries:
        spec = importlib.util.find_spec(lib)
        if spec is None:
            raise ImportError(
                f"Missing Python library '{lib}'. Please install it, e.g. pip install {lib}."
            )


def which_executable(exec_name, tool_label=None):
    """
    Returns the absolute path to 'exec_name' if found in PATH or None otherwise.
    If 'tool_label' is given, it is used in the error message to better identify 
    the missing tool. 
    """
    path = shutil.which(exec_name)
    if path is None and tool_label:
        raise FileNotFoundError(
            f"Missing required tool '{tool_label}' (executable: '{exec_name}'). "
            f"Install or provide a path so that '{exec_name}' is in your $PATH."
        )
    return path


def check_dependencies(args):
    """
    Inspect args to see which external tools are *actually needed*, 
    then verify that they are available. 
    Raise an exception with a useful error message if something is missing.
    """

    logger = logging.getLogger("galba_pipeline.deps")
    
    # 1) Always required: diamond, miniprot, miniprot_boundary_scorer, TSEBRA, AUGUSTUS, etc.
    #    We'll do a quick check that these are in PATH or a user-specified location:
    #    If the script eventually calls them via a path like args.DIAMOND_PATH, 
    #    you might check that path. For simplicity here, we check in $PATH:

    # Example approach:
    which_executable("diamond", "DIAMOND")
    which_executable("miniprot", "miniprot")
    which_executable("miniprot_boundary_scorer", "miniprot_boundary_scorer")
    which_executable("tsebra.py", "TSEBRA")      # Or just "tsebra" if that’s how your pipeline calls it.
    
    # For AUGUSTUS, you might check for 'augustus' and 'etraining'.
    # (You also might have a dedicated config path. This example simply 
    #  checks whether it’s in PATH.)
    which_executable("augustus", "AUGUSTUS (augustus)")
    which_executable("etraining", "AUGUSTUS (etraining)")
    which_executable("optimize_augustus.pl", "optimize_augustus.pl")
    which_executable("new_species.pl", "new_species.pl")
    
    # 2) Python libraries
    check_python_libraries()

    # 3) Conditionals for RNA-seq / Iso-Seq 
    #    If user provided FASTQ-based short reads (rnaseq_fq), we need HISAT2:
    if args.rnaseq_fq and len(args.rnaseq_fq) > 0:
        logger.info("User provided RNA-seq in FASTQ => Checking for 'hisat2' tools.")
        which_executable("hisat2", "HISAT2")
        which_executable("hisat2-build", "HISAT2 (hisat2-build)")
    
    # If user provided FASTQ-based Iso-Seq reads => minimap2 is required:
    if args.isoseq_fq and len(args.isoseq_fq) > 0:
        logger.info("User provided Iso-Seq in FASTQ => Checking for 'minimap2'.")
        which_executable("minimap2", "minimap2")
    
    # If user provided BAM-based RNA-seq or Iso-Seq, we still need samtools unless user gave a StringTie file
    # but from your statement, you said "for any rnaseq or isoseq that is not stringtie, samtools is needed."
    # So if we see either isoseq_bam or rnaseq_bam, or rnaseq_fq or isoseq_fq, we need samtools:
    need_samtools = False
    if (args.rnaseq_bam and len(args.rnaseq_bam) > 0) or \
       (args.isoseq_bam and len(args.isoseq_bam) > 0) or \
       (args.rnaseq_fq and len(args.rnaseq_fq) > 0) or \
       (args.isoseq_fq and len(args.isoseq_fq) > 0):
        
        if not args.stringtie_gtf:
            need_samtools = False
    
    if need_samtools:
        logger.info("RNA-seq/Iso-Seq data provided that is not StringTie => Checking for 'samtools'.")
        which_executable("samtools", "samtools")

    # 4) For any RNASEQ or ISOSEQ that is not a *pre-assembled* stringtie,
    #    we also need stringtie to do the transcript assembly:
    need_stringtie = False
    if not args.stringtie_gtf:
        # means user must do assembly themselves => stringtie is needed
        if (args.rnaseq_fq and len(args.rnaseq_fq) > 0) or \
           (args.rnaseq_bam and len(args.rnaseq_bam) > 0) or \
           (args.isoseq_fq and len(args.isoseq_fq) > 0) or \
           (args.isoseq_bam and len(args.isoseq_bam) > 0):
            need_stringtie = True
    
    if need_stringtie:
        logger.info("No StringTie GTF => Need to run StringTie ourselves.")
        which_executable("stringtie", "StringTie")

    # 5) If any rnaseq or isoseq was provided, we need TransDecoder + bedtools
    if (args.rnaseq_fq or args.rnaseq_bam or 
        args.isoseq_fq or args.isoseq_bam or 
        args.stringtie_gtf):
        logger.info("Transcript data was provided => Checking TransDecoder + bedtools.")
        which_executable("TransDecoder.LongOrfs", "TransDecoder.LongOrfs")
        which_executable("TransDecoder.Predict", "TransDecoder.Predict")
        # We also call the TransDecoder "util" scripts, so either ensure they are in PATH
        # or do additional checks if your pipeline runs them by name:
        # e.g. check "gtf_genome_to_cdna_fasta.pl" if needed
        which_executable("bedtools", "bedtools")  # Usually just "bedtools".
    
    logger.info("All required tools appear to be available for the given pipeline configuration.")
