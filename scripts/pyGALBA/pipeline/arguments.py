# pipeline/arguments.py
import argparse

def parse_arguments(cli_args):
    parser = argparse.ArgumentParser(
        description="GALBA is a fully automated pipeline for the annotation of protein coding genes."
    )
    
    # Existing arguments
    parser.add_argument("--genome", required=True, help="Path to genome FASTA file.")
    parser.add_argument("--prot_seq", help="Protein FASTA file for alignment.")
    parser.add_argument("--hints", nargs='*', default=[],
                        help="Additional GFF hints file(s).")
    parser.add_argument("--species", help="Species name for AUGUSTUS.")
    parser.add_argument("--workingdir", default=".",
                        help="Path to working directory for output files.")
    parser.add_argument("--crf", action="store_true", help="Enable CRF training.")
    parser.add_argument("--nice", action="store_true", help="Run external tools with 'nice'.")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use.")
    
    # The new input options
    parser.add_argument("--isoseq_fq", nargs='*',
        help="One or more FASTQ files containing Iso-Seq reads.")
    parser.add_argument("--rnaseq_fq", nargs='*',
        help="One or more FASTQ files containing RNA-seq reads.")
    parser.add_argument("--isoseq_bam", nargs='*',
        help="One or more BAM files containing Iso-Seq alignments.")
    parser.add_argument("--rnaseq_bam", nargs='*',
        help="One or more BAM files containing RNA-seq alignments.")
    parser.add_argument("--stringtie_gtf",
        help="StringTie GTF output providing transcripts for evidence.")
    
    # Possibly additional flags...
    
    # Additional advanced or debugging flags
    parser.add_argument("--skipAllTraining", action="store_true")
    parser.add_argument("--skipOptimize", action="store_true")
    
    args = parser.parse_args(cli_args)
    return args
