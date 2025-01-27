# pipeline/fileformats.py

import os
import sys
from Bio import SeqIO

def check_fileformats(args):
    """
    Check that each input file is present and, if the pipeline requires 
    a specific format, attempt minimal format verification.
    
    Accumulate all format or file problems in a list of error messages, 
    then raise a single exception if needed.
    """

    errors = []

    # 1) Genome (FASTA)
    if not os.path.isfile(args.genome):
        errors.append(f"Genome file not found: {args.genome}")
    else:
        # Minimal check for FASTA
        if not is_fasta(args.genome):
            errors.append(
                f"Genome file {args.genome} does not appear to be a FASTA file."
            )
        # Check if it's DNA or protein
        seq_type = dna_or_protein(args.genome)
        if seq_type == "Protein":
            errors.append(
                f"Genome file {args.genome} appears to contain protein sequences."
            )

    # 2) Protein sequences file (FASTA)
    if not os.path.isfile(args.prot_seq):
        errors.append(f"Protein file not found: {args.prot_seq}")
    else:
        if not is_fasta(args.prot_seq):
            errors.append(
                f"Protein file {args.prot_seq} does not appear to be a FASTA file."
            )
        # Check if it's DNA or protein
        seq_type = dna_or_protein(args.prot_seq)
        if seq_type == "DNA":
            errors.append(
                f"Protein file {args.prot_seq} appears to contain DNA sequences."
            )
    
    # 3) If user gave any rnaseq_fq => check each is present & is FASTQ
    if args.rnaseq_single_fq or args.rnaseq_paired_fq:
        for fq in args.rnaseq_single_fq + args.rnaseq_paired_fq:
            if not os.path.isfile(fq):
                errors.append(f"RNA-seq FASTQ file not found: {fq}")
            else:
                if not is_fastq(fq):
                    errors.append(
                        f"RNA-seq file {fq} does not appear to be FASTQ."
                    )

    # 4) If user gave any isoseq_fq => check each is present & is FASTQ
    if args.isoseq_fq:
        for fq in args.isoseq_fq:
            if not os.path.isfile(fq):
                errors.append(f"Iso-Seq FASTQ file not found: {fq}")
            else:
                if not is_fastq(fq):
                    errors.append(
                        f"Iso-Seq file {fq} does not appear to be FASTQ."
                    )
    
    # 5) If user gave any rnaseq_bam => check each is present & is BAM
    if args.rnaseq_bam:
        for bam in args.rnaseq_bam:
            if not os.path.isfile(bam):
                errors.append(f"RNA-seq BAM file not found: {bam}")
            else:
                if not is_bam(bam):
                    errors.append(
                        f"RNA-seq file {bam} does not appear to be BAM format."
                    )

    # 6) If user gave any isoseq_bam => check each is present & is BAM
    if args.isoseq_bam:
        for bam in args.isoseq_bam:
            if not os.path.isfile(bam):
                errors.append(f"Iso-Seq BAM file not found: {bam}")
            else:
                if not is_bam(bam):
                    errors.append(
                        f"Iso-Seq file {bam} does not appear to be BAM format."
                    )

    # 7) If user gave a stringtie_gtf => check presence & minimal GTF check
    if args.stringtie_gtf:
        if not os.path.isfile(args.stringtie_gtf):
            errors.append(f"StringTie GTF file not found: {args.stringtie_gtf}")
        else:
            if not is_gtf(args.stringtie_gtf):
                errors.append(
                    f"StringTie file {args.stringtie_gtf} does not appear to be GTF."
                )

    # If any errors accumulated, raise a single exception
    if errors:
        err_msg = "One or more file format issues were detected:\n  - "
        err_msg += "\n  - ".join(errors)
        sys.exit(err_msg)

    # Otherwise, do nothing
    # Print info if desired
    # print("All file formats look OK.")

def is_fasta(filepath):
    """
    Minimal check for FASTA: read first line(s) => starts with '>'
    Returns True if likely FASTA, else False.
    """
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
        return first_line.startswith(">")
    except:
        return False

def is_fastq(filepath):
    """
    Minimal check for FASTQ: 
    read first few lines => line1 starts with '@', line3 starts with '+'
    """
    try:
        with open(filepath, 'r') as f:
            lines = []
            for _ in range(4):
                line = f.readline()
                if not line:
                    break
                lines.append(line.strip())
        if len(lines) < 3:
            return False
        # check line 1 starts with '@' and line 3 starts with '+'
        if lines[0].startswith('@') and lines[2].startswith('+'):
            return True
        return False
    except:
        return False

def is_bam(filepath):
    """
    Minimal check for BAM. 
    Typically you'd do: read first 4 bytes => check for 'BAM\1'
    But for simplicity, we do a naive approach or extension check.
    """
    # If extension isn't .bam => fail
    # A more robust approach is to open the file in binary mode 
    # and check for 'BAM\1'
    if not filepath.lower().endswith(".bam"):
        return False
    # Optionally read the magic bytes
    try:
        with open(filepath, "rb") as f:
            magic = f.read(4)
        return magic == b'BAM\1'
    except:
        return False

def is_gtf(filepath):
    """
    Minimal check for GTF. 
    We check if extension is .gtf or .gff, 
    then check the first line for tab-delimited & a possible 'transcript_id'.
    """
    extension = os.path.splitext(filepath)[1].lower()
    if extension not in [".gtf", ".gff", ".gff3"]:
        return False
    
    # read a couple lines to see if it has e.g. 9 fields, 
    # or a 'transcript_id' in the 9th
    try:
        with open(filepath, 'r') as f:
            for _ in range(10):
                line = f.readline()
                if not line:
                    break
                if line.startswith("#"):
                    continue
                parts = line.strip().split('\t')
                if len(parts) == 9:
                    # Possibly check if 'transcript_id' in parts[-1]
                    if 'transcript_id' in parts[-1]:
                        return True
            return False
    except:
        return False

def dna_or_protein(filepath):
    """
    Check if a fasta file is DNA or protein sequences.
    We'll check the first line(s) to see if it's DNA or protein.
    Use BioPython for a robust check.
    """
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
        if first_line.startswith(">"):
            # Use BioPython to parse the file and check the first sequence
            seq = next(SeqIO.parse(filepath, "fasta"))
            # Check if it's a DNA or protein sequence
            if set(str(seq.seq).upper()) <= set("ACGTN"):
                return "DNA"
            return "Protein"
    except:
        return "Unknown"