#!/usr/bin/env python3
"""
Author: Katharina J. Hoff and Mario Stanke
Ported to Python3 by: Katharina J. Hoff

This script samples training gene structures according to their intron number
(originally from GeneMark-EX). Downsampling single-exon genes and genes with
fewer introns (1–5) can improve the accuracy of AUGUSTUS when trained on these
downsampled structures.

Equivalent of:
  downsample_traingenes.pl --in_gtf=traingenes.gtf --out_gtf=out.gtf [options]

Options:
    --in_gtf            Input GTF file of training gene structures
    --out_gtf           Output GTF file with downsampled structures
    --lambda            Poisson parameter lambda (default 0)
    --intron_num_lst    Output a file listing <intron_count> <transcript_id>
    --version           Print the script version (1.0)
    --help              Print this help message

Example:
    python3 downsample_traingenes.py --in_gtf=traingenes.gtf --out_gtf=out.gtf --lambda=1
"""
import sys
import re
import math
import random
import argparse

__version__ = "1.0.0"

def poisson_pmf(k: int, lam: float) -> float:
    """Compute the Poisson probability mass function P(X = k) with parameter lam."""
    # P(X=k) = e^-lam * (lam^k / k!)
    return (math.e ** (-lam)) * (lam ** k) / math.factorial(k)

def main():
    parser = argparse.ArgumentParser(
        description="Downsample training gene structures from GeneMark-ET/EP/ETP."
    )
    parser.add_argument(
        "--in_gtf", required=True, help="Input GTF file with training gene structures."
    )
    parser.add_argument(
        "--out_gtf", required=True, help="Output GTF file with downsampled gene structures."
    )
    parser.add_argument(
        "--lambda", type=float, default=0,
        help="Parameter lambda of the Poisson distribution (default: 0)."
    )
    parser.add_argument(
        "--intron_num_lst",
        help="Optional: file to write <intron_count> <transcript_id> for each kept gene."
    )
    parser.add_argument(
        "--version", action="store_true", help="Print version of the script and exit."
    )

    args = parser.parse_args()

    # If --version is requested, print and exit.
    if args.version:
        print(f"Version {__version__}")
        sys.exit(0)

    in_gtf = args.in_gtf
    out_gtf = args.out_gtf
    lam = args.__dict__["lambda"]  # retrieve the value of --lambda
    intron_num_lst = args.intron_num_lst

    # Make downsampling reproducible (similar to 'srand 10' in Perl)
    random.seed(10)

    # Data structures to store transcripts
    nIntrons = {}   # transcript_id -> (count of introns)
    tx_lines = {}   # transcript_id -> list of lines (GTF entries)

    # Read the GTF file and parse introns/CDS
    # The original Perl logic: each time a line has 'CDS', we do nIntrons[txid]++.
    # So intron count = #CDS lines - 1 for each transcript, effectively.
    with open(in_gtf, "r") as f_in:
        for line in f_in:
            match = re.search(r'transcript_id\s+"(\S+)"', line)
            if match:
                txid = match.group(1)
                if txid not in tx_lines:
                    tx_lines[txid] = []
                tx_lines[txid].append(line)

                # If line is a CDS, increment the intron count (just like in the Perl script)
                if "\tCDS\t" in line:
                    if txid not in nIntrons:
                        nIntrons[txid] = 0
                    else:
                        nIntrons[txid] += 1

    # Compute the cumulative distribution function (CDF) F of Poisson(λ) for i=0..5
    max_intron_number = 5
    F = []
    for i in range(max_intron_number + 1):
        if i == 0:
            F.append(poisson_pmf(i, lam))
        else:
            F.append(F[i-1] + poisson_pmf(i, lam))

    # Prepare outputs
    out_handle = open(out_gtf, "w")
    lst_handle = None
    if intron_num_lst:
        lst_handle = open(intron_num_lst, "w")

    # Keep track of single-exon genes, ensure at least 20 are kept
    min_single_exon_genes = 20
    single_exon_gene_counter = 0

    # For not writing the same intronNum multiple times in intron_num_lst
    intronNumPrinted = set()

    # Sort transcript IDs for stable output (like the Perl script)
    for txid in sorted(nIntrons.keys()):
        intronNum = nIntrons[txid]
        if intronNum == 0:
            single_exon_gene_counter += 1

        # Draw a uniform random number
        u = random.random()
        # If intronNum > max_intron_number, we always keep the gene
        # OR if u <= F[intronNum], we keep it
        # OR if intronNum==0 and we haven't kept 20 single-exon yet
        # then we keep the gene. This replicates the original sampling logic.
        if (intronNum > max_intron_number
                or u <= F[intronNum]
                or (intronNum == 0 and single_exon_gene_counter <= min_single_exon_genes)):
            for gtf_line in tx_lines[txid]:
                out_handle.write(gtf_line)
            if lst_handle and txid not in intronNumPrinted:
                lst_handle.write(f"{intronNum}\t{txid}\n")
                intronNumPrinted.add(txid)

    out_handle.close()
    if lst_handle:
        lst_handle.close()

    # Print a warning if fewer than 20 single-exon genes are found overall
    if single_exon_gene_counter < 20:
        print(f"WARNING: Number of single-exon training genes is smaller than 20. "
              f"It is: {single_exon_gene_counter}!",
              file=sys.stdout)

if __name__ == "__main__":
    main()
