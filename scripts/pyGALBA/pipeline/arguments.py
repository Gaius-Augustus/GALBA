# pipeline/arguments.py
import argparse

def parse_arguments(cli_args):
    parser = argparse.ArgumentParser(
        description="GALBA is a fully automated pipeline for the annotation of protein coding genes."
    )
    
    # Always required:
    parser.add_argument("--genome", required=True,
        help="Path to genome FASTA file.")
    parser.add_argument("--prot_seq", required=True,
        help="Protein FASTA file for alignment (e.g. reference proteome).")
    
    # Additional hints:
    parser.add_argument("--hints", nargs='*', default=[],
        help="Additional GFF hints file(s).")
    
    # Standard arguments:
    parser.add_argument("--species",
        help="Species name for AUGUSTUS.")
    parser.add_argument("--workingdir", default="./GALBA_output",
        help="Path to working directory for output files.")
    parser.add_argument("--crf", action="store_true",
        help="Enable CRF training.")
    parser.add_argument("--nice", action="store_true",
        help="Run external tools with 'nice'.")
    parser.add_argument("--threads", type=int, default=1,
        help="Number of threads to use.")
    
    # New RNA-seq / Iso-Seq input options; all optional, but restricted:
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
    
    # Protein database required if any RNA-Seq/Iso-Seq/StringTie evidence given
    parser.add_argument("--protein_db",
        help="Large protein database in FASTA format (not the same as --prot_seq). "
             "Required if any Iso-Seq/RNA-seq/StringTie input is provided.")
    
    # Additional advanced or debugging flags
    parser.add_argument("--skipAllTraining", action="store_true",
        help="Skip all AUGUSTUS training steps.")
    parser.add_argument("--skipOptimize", action="store_true",
        help="Skip parameter optimization steps for AUGUSTUS.")
    
    args = parser.parse_args(cli_args)

    # -------------------------------------------------------------------------
    # Validate the RNA-seq / Iso-Seq input logic:
    # The user can select EXACTLY one 'mode' of using transcript evidence:
    #   1) FASTQ files (isoseq_fq and/or rnaseq_fq),
    #   2) BAM files (isoseq_bam and/or rnaseq_bam),
    #   3) StringTie GTF,
    #   4) Or no transcript evidence at all (use no iso/rna/stringtie).
    # -------------------------------------------------------------------------

    has_fastq = bool(args.isoseq_fq or args.rnaseq_fq)
    has_bam = bool(args.isoseq_bam or args.rnaseq_bam)
    has_stringtie = bool(args.stringtie_gtf)

    modes_selected = sum([has_fastq, has_bam, has_stringtie])
    if modes_selected > 1:
        parser.error(
            "Invalid combination of RNA-seq / Iso-Seq inputs. "
            "Only one mode is allowed: either (fastq) or (bam) or (stringtie)."
        )
    
    # If the user has chosen ANY transcript evidence, then --protein_db is required
    if modes_selected == 1:
        if not args.protein_db:
            parser.error(
                "You must supply --protein_db whenever you provide "
                "Iso-Seq / RNA-seq / StringTie evidence."
            )

    return args
