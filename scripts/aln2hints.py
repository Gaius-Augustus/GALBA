#!/usr/bin/env python3
"""
aln2hints.py

This script generates hints from alignment outputs (miniprot or GenomeThreader (gth)).

Ported from the original Perl script aln2hints.pl by Katharina Hoff (November 4, 2022).

Usage (examples):
  aln2hints.py --in=align.gff --out=hintsfile.gff --prg=miniprot
  aln2hints.py --in=align.gff --out=hintsfile.gff --prg=gth
  aln2hints.py --help

Arguments:
  --in                Input alignment file (gth: GFF3, miniprot: GTF).
  --out               Output file containing hints (CDSpart/CDS, intron, start, stop).
  --prg               Alignment program: 'gth' or 'miniprot'.

Options:
  --CDSpart_cutoff=n  Number of base pairs to cut off each CDSpart hint (default 15).
  --maxintronlen=n    Maximum allowed intron length (default 350000).
  --minintronlen=n    Minimum allowed intron length (default 41).
  --priority=n        Priority of the generated hints (default 4).
  --CDS               If set, output complete CDS hints instead of partial (CDSpart).
  --genome_file=s     If prg=miniprot and you want start/stop hints, specify the genome FASTA
                      from which the alignments were generated.
  --source=s          Source identifier for all hints (default 'P').
  --dir=s             A working directory path (mainly for internal usage).
  --version           Print script version.
  --help              Print this help message.

Example:
  aln2hints.py --in=align.gtf --out=hintsfile.gff --prg=miniprot --genome_file=genome.fa
"""

import sys
import os
import re
import math
import argparse

VERSION = "1.0"


def read_genome_fasta(genome_file):
    """Read a FASTA file into a dict: {header_without_gt: sequence_string}."""
    genome_dict = {}
    current_header = None

    with open(genome_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_header = line[1:]  # remove '>'
                genome_dict[current_header] = ""
            else:
                if current_header is not None:
                    genome_dict[current_header] += line
    return genome_dict


def main():
    parser = argparse.ArgumentParser(
        add_help=False,
        description=("Generate hints from miniprot or GenomeThreader (gth) output.")
    )
    parser.add_argument("--help", action="store_true", help="Print help message.")
    parser.add_argument("--version", action="store_true", help="Print script version.")
    parser.add_argument("--in", dest="alignfile", help="Alignment input file.")
    parser.add_argument("--out", dest="hintsfilename", help="Output hints file.")
    parser.add_argument("--prg", dest="prgsrc", choices=["gth", "miniprot"],
                        help="Alignment program: 'gth' or 'miniprot'.")
    parser.add_argument("--CDSpart_cutoff", type=int, default=15,
                        help="Number of base pairs to cut off each CDSpart (default 15).")
    parser.add_argument("--maxintronlen", type=int, default=350000,
                        help="Max allowed intron length (default 350000).")
    parser.add_argument("--minintronlen", type=int, default=41,
                        help="Min allowed intron length (default 41).")
    parser.add_argument("--priority", type=int, default=4,
                        help="Priority for all generated hints (default 4).")
    parser.add_argument("--source", default="P",
                        help="Source identifier (default 'P').")
    parser.add_argument("--dir", help="Working directory path.")
    parser.add_argument("--CDS", action="store_true",
                        help="If set, output complete CDS hints instead of partial (CDSpart).")
    parser.add_argument("--genome_file",
                        help="Genome FASTA for generating start/stop hints (used with miniprot).")

    args = parser.parse_args()

    # Handle --help
    if args.help:
        print(__doc__)
        parser.print_help()
        sys.exit(0)

    # Handle --version
    if args.version:
        print(f"aln2hints.py version {VERSION}")
        sys.exit(0)

    # Basic checks
    if not args.alignfile or not args.hintsfilename or not args.prgsrc:
        sys.stderr.write(
            "ERROR: Missing required arguments --in, --out, or --prg.\n\n"
        )
        print(__doc__)
        parser.print_help()
        sys.exit(1)

    # Possibly set dir
    directory = args.dir
    if not directory:
        directory = os.getcwd()
    else:
        directory = os.path.abspath(directory)
    if directory.endswith("/"):
        directory = directory[:-1]

    # Check existence of input file
    alignfile = os.path.abspath(args.alignfile)
    if not os.path.exists(alignfile):
        sys.stderr.write(
            f"ERROR: Alignment file {alignfile} does not exist. Please check.\n"
        )
        sys.exit(1)

    # Determine prgsrc: 'gth2h' or 'miniprot2h'
    prgsrc = None
    intron_threshold = None
    intron_threshold_gth = 0.7
    if args.prgsrc == "gth":
        prgsrc = "gth2h"
        intron_threshold = intron_threshold_gth
    elif args.prgsrc == "miniprot":
        prgsrc = "miniprot2h"

    # If genome file is provided but prg != miniprot2h, ignore
    genome_dict = {}
    if prgsrc == "miniprot2h" and args.genome_file:
        if not os.path.exists(args.genome_file):
            sys.stderr.write(f"ERROR: genome_file {args.genome_file} does not exist.\n")
            sys.exit(1)
        # Read genome if miniprot2h with a provided genome_file
        genome_dict = read_genome_fasta(args.genome_file)
    elif prgsrc != "miniprot2h" and args.genome_file:
        sys.stderr.write(
            f"ERROR: program name is {prgsrc} but a genome_file was specified. "
            "Ignoring genome file.\n"
        )

    # Decide if we use "CDS" or "CDSpart"
    if args.CDS:
        CDSpartid = "CDS"
    else:
        CDSpartid = "CDSpart"

    # Keep track of variables equivalent to the Perl script's
    parent = None
    prevParent = "noP"
    prevScore = 0.0
    intron_start = None
    intron_end = None

    # Hard-coded for scipio logic in original, not used but we keep placeholders
    qstart = None
    qend = None
    prevQend = None

    # open input, output
    try:
        aln_in = open(alignfile, "r")
    except OSError:
        sys.stderr.write(f"Cannot open file: {alignfile}\n")
        sys.exit(1)

    hints_outfile = os.path.abspath(args.hintsfilename)
    try:
        hints_out = open(hints_outfile, "w")
    except OSError:
        sys.stderr.write(f"Cannot open file: {hints_outfile}\n")
        sys.exit(1)

    # Helper functions nested here for direct access to local variables
    def print_hint_line(seqname, src, ftype, start_, end_, score_, strand_, frame_="."):
        """Write a single hint line to output in GFF-like format."""
        line = (
            f"{seqname}\t{src}\t{ftype}\t{start_}\t{end_}\t{score_}\t"
            f"{strand_}\t{frame_}\tsrc={args.source};grp={parent};pri={args.priority}\n"
        )
        hints_out.write(line)

    def print_start(seqname, strand, start_, end_, prnt):
        """Replicate 'print_start' subroutine from Perl."""
        # start codon is 3 bp
        # If forward: start..start+2
        # If reverse: end-2..end
        ftype = "start"
        if strand == "+":
            s_ = start_
            e_ = start_ + 2
        else:
            s_ = end_ - 2
            e_ = end_
        line = (
            f"{seqname}\t{prgsrc}\t{ftype}\t{s_}\t{e_}\t.\t{strand}\t0\t"
            f"src={args.source};grp={prnt};pri={args.priority}\n"
        )
        hints_out.write(line)

    def print_stop(seqname, strand, start_, end_):
        """Replicate 'print_stop' subroutine from Perl."""
        # stop codon is 3 bp
        # If forward: end-2..end
        # If reverse: start..start+2
        ftype = "stop"
        if strand == "+":
            s_ = end_ - 2
            e_ = end_
        else:
            s_ = start_
            e_ = start_ + 2
        line = (
            f"{seqname}\t{prgsrc}\t{ftype}\t{s_}\t{e_}\t.\t{strand}\t0\t"
            f"src={args.source};grp={parent};pri={args.priority}\n"
        )
        hints_out.write(line)

    def write_intron_line(f_array):
        """Equivalent to 'get_intron' subroutine in Perl."""
        nonlocal intron_start, intron_end, prevScore, prevParent, qstart, qend, prevQend
        nonlocal parent

        line_score = 0.0
        # gth2h logic uses an average of previous + current
        if prgsrc == "gth2h":
            line_score = prevScore + float(f_array[5]) / 2.0
        else:
            line_score = 0.0

        strand_ = f_array[6]
        seqname_ = f_array[0]

        # If we have a new parent, we set up the new intron boundary
        if prevParent != parent:
            if strand_ == "-" and prgsrc == "miniprot2h":
                # aligners output in reverse order
                # so intron_end is the next region's start - 1
                intron_end = int(f_array[3]) - 1
            else:
                intron_start = int(f_array[4]) + 1
        else:
            # We have the second exon/CDS line in a row => we can define an intron
            if strand_ == "-" and prgsrc == "miniprot2h":
                intron_start = int(f_array[4]) + 1
            else:
                intron_end = int(f_array[3]) - 1

            # Ensure intron_start <= intron_end
            if intron_end < intron_start:
                temp = intron_start
                intron_start = intron_end
                intron_end = temp

            if prgsrc == "miniprot2h":
                # Print intron line
                print_hint_line(
                    seqname_, prgsrc, "intron",
                    intron_start, intron_end,
                    line_score, strand_, "."
                )
            else:
                # For gth, check constraints: length >= minintronlen, <= maxintronlen,
                # and intron_score > intron_threshold if threshold is defined
                length_ = intron_end - intron_start + 1
                if (length_ >= args.minintronlen and
                        length_ <= args.maxintronlen and
                        (intron_threshold is None or line_score > intron_threshold)):
                    # Additional gth condition: scipio logic (not used here, but we replicate)
                    # ( defined($prevQend) && $prevQend + 1 == $qstart ) in original
                    # No actual qstart usage in the miniprot/gth code in the original,
                    # so we mimic the check if needed. Here we'll skip the scipio check.
                    print_hint_line(
                        seqname_, prgsrc, "intron",
                        intron_start, intron_end,
                        line_score, strand_, "."
                    )

            # Prepare for next intron boundary
            if strand_ == "-" and prgsrc == "miniprot2h":
                intron_end = int(f_array[3]) - 1
            else:
                intron_start = int(f_array[4]) + 1

        prevScore = float(f_array[5]) / 2.0
        prevParent = parent

    # Main loop over lines
    for line in aln_in:
        if "*" in line:
            # skip lines containing a star
            continue
        line = line.strip()
        f = line.split("\t")
        if len(f) < 8:
            continue

        seqname = f[0]
        ftype = f[2]
        start_ = int(f[3])
        end_ = int(f[4])
        score_ = f[5]
        strand_ = f[6]
        frame_ = f[7]
        attr = f[8] if len(f) > 8 else ""

        # Make sure start <= end
        if end_ < start_:
            start_, end_ = end_, start_

        # For gth2h, parent is extracted from lines with type == "mRNA"
        if ftype == "mRNA" and prgsrc == "gth2h":
            # example: seq=chr1; ID=rna1; geneID=gene1 ...
            # The Perl version does something like:
            #   my @info   = split(/\=/, $f[8]);
            #   my @rnaid  = split(/;/, $info[1]);
            #   ...
            # Then constructs parent = info[0]_<seqname>_<rnaid[0]>_<geneid[0]>.
            # We mimic that logic approximately:
            # But gth often has lines like:
            #     .   .   mRNA    .   .   .   .   .   ID=mRNA1;rnaID=xxx;geneID=xxx
            # We'll do minimal parsing consistent with the original script.
            info_parts = attr.split("=")
            # This is a heuristic approach to reconstruct the parent:
            if len(info_parts) > 2:
                # e.g. "ID=rna1; geneID=gene1"
                rnaid_val = info_parts[1].split(";")[0]
                geneid_val = info_parts[2].split(";")[0]
                # final piece might have trailing spaces:
                rnaid_val = rnaid_val.strip()
                geneid_val = geneid_val.strip()
                # Rebuild parent
                parent = f"{info_parts[0]}_{seqname}_{rnaid_val}_{geneid_val}"
            else:
                # fallback if unexpected
                parent = attr
            continue

        # For miniprot2h, parent is extracted from lines with type == "gene"
        if ftype == "gene" and prgsrc == "miniprot2h":
            # The Perl script does something like: 
            # f[8] =~ m/gene_id "([^"]+)"/ => parent = group(1)
            match = re.search(r'gene_id\s+"([^"]+)"', attr)
            if match:
                parent = match.group(1)
            else:
                parent = "unknownParent"

            # If genome_file is provided => create start & stop hints
            if args.genome_file:
                # Start codon check
                pot_start_seq = None
                if strand_ == "+":
                    # substring: start_-1 .. start_-1+2
                    if seqname in genome_dict:
                        # watch out for indexing issues
                        # Python substring: [start_-1:start_-1+3]
                        pot_start_seq = genome_dict[seqname][(start_ - 1): (start_ - 1 + 3)]
                else:
                    # substring: end_-3 .. end_-1 reversed/complemented
                    if seqname in genome_dict:
                        pot_start_seq = genome_dict[seqname][(end_ - 3):end_]
                        # complement, then reverse
                        # a->t, c->g, g->c, t->a
                        pot_start_seq = pot_start_seq.translate(
                            str.maketrans("acgtACGT", "tgcaTGCA")
                        )
                        pot_start_seq = pot_start_seq[::-1]

                if pot_start_seq:
                    if re.match(r'(ATG|TTG|GTG|CTG)', pot_start_seq, re.IGNORECASE):
                        print_start(seqname, strand_, start_, end_, parent)

                # Stop codon check
                pot_stop_seq = None
                if strand_ == "+":
                    # substring: end_ .. end_+2
                    if seqname in genome_dict and end_ + 2 < len(genome_dict[seqname]):
                        pot_stop_seq = genome_dict[seqname][end_: end_ + 3]
                else:
                    # substring: start_-4 .. start_-2 reversed/complemented
                    if seqname in genome_dict and (start_ - 4) >= 0:
                        pot_stop_seq = genome_dict[seqname][(start_ - 4):(start_ - 1)]
                        pot_stop_seq = pot_stop_seq.translate(
                            str.maketrans("acgtACGT", "tgcaTGCA")
                        )
                        pot_stop_seq = pot_stop_seq[::-1]

                if pot_stop_seq:
                    if re.match(r'(TAA|TGA|TAG)', pot_stop_seq, re.IGNORECASE):
                        print_stop(seqname, strand_, start_, end_)

            continue

        # For lines with type == "CDS" or "cds" or "protein_match":
        if ftype in ["CDS", "cds", "protein_match"]:
            # If we do not output the entire CDS, we cut off each end by CDSpart_cutoff
            start_adj = start_
            end_adj = end_
            if not args.CDS:
                start_adj = start_ + args.CDSpart_cutoff
                end_adj = end_ - args.CDSpart_cutoff
                if start_adj > end_adj:
                    # center them (like the Perl code)
                    midpoint = (start_adj + end_adj) // 2
                    start_adj = end_adj = midpoint

            print_hint_line(
                seqname, prgsrc, CDSpartid,
                start_adj, end_adj,
                score_, strand_, frame_
            )

            # For miniprot2h, also attempt intron logic
            if prgsrc == "miniprot2h":
                write_intron_line(f)

        # For lines with type == "exon" and prgsrc == "gth2h":
        if ftype == "exon" and prgsrc == "gth2h":
            write_intron_line(f)

    aln_in.close()
    hints_out.close()


if __name__ == "__main__":
    main()
