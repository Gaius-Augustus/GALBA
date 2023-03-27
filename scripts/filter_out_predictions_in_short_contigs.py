#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(
    description='Delete predictions in short contigs from gtf output file for GALBA (or BRAKER)')
parser.add_argument('-f', '--gtf_file', type=str, required=True,
                    help="Input gtf file (from GALBA or BRAKER)")
parser.add_argument('-g', '--genome_file', type=str, required=True,
                    help="Input genome file in fasta format")
parser.add_argument('-l', '--min_length', type=int, default=5000,
                    help="Minimum length of contigs to keep predictions in")
args = parser.parse_args()

# read genome file into a dictionary that has contig name as key and contig length as value
genome_dict = {}
try:
    with open(args.genome_file, "r") as genome_file:
        for line in genome_file:
            # remove trailing newline
            line = line.rstrip()
            if re.search(r'^>', line):
                header = re.search(r'^>([^ ]+)', line).group(1)
                genome_dict[header] = 0
            else:
                genome_dict[header] += len(line.strip())
except IOError:
    print("Could not read file:", args.genome_file)
    sys.exit(1)

# loop through gtf file and print lines that are not in short contigs
try:
    with open(args.gtf_file, "r") as gtf_file:
        for line in gtf_file:
            if genome_dict[re.split(r'\t', line)[0]] >= args.min_length:
                print(line.strip())
except IOError:
    print("Could not read file:", args.gtf_file)
    sys.exit(1)
