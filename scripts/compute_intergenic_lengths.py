#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(
    description='Compute intergenic regions between gene features in gtf file')
parser.add_argument('-f', '--gtf_file', type=str, required=True,
                    help="Input gtf file (from GALBA or BRAKER)")
args = parser.parse_args()


def compute_intergenic_lengths(gtf_file):
    gene_positions = {}  # dictionary to store gene positions
    intergenic_lengths = []  # list to store intergenic lengths

    # Parse the GTF file
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            feature_type = parts[2]
            if feature_type != 'gene':
                continue

            seqname = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            # Store gene positions based on the sequence name and strand
            if seqname not in gene_positions:
                gene_positions[seqname] = {}

            if strand not in gene_positions[seqname]:
                gene_positions[seqname][strand] = []

            gene_positions[seqname][strand].append((start, end))

    # Compute intergenic lengths
    for seqname in gene_positions:
        if '+' in gene_positions[seqname] and '-' in gene_positions[seqname]:
            # Sort gene positions based on the start coordinate
            gene_positions[seqname]['+'].sort(key=lambda x: x[0])
            gene_positions[seqname]['-'].sort(key=lambda x: x[0])

            # Compute intergenic lengths between adjacent genes on the same strand
            for i in range(len(gene_positions[seqname]['+']) - 1):
                intergenic_length = gene_positions[seqname]['+'][i + 1][0] - gene_positions[seqname]['+'][i][1] - 1
                intergenic_lengths.append(intergenic_length)

            for i in range(len(gene_positions[seqname]['-']) - 1):
                intergenic_length = gene_positions[seqname]['-'][i + 1][0] - gene_positions[seqname]['-'][i][1] - 1
                intergenic_lengths.append(intergenic_length)

    return intergenic_lengths


intergenic_lengths = compute_intergenic_lengths(args.gtf_file)
print("Intergenic lengths:", (sum(intergenic_lengths)/len(intergenic_lengths)))
