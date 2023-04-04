#!/usr/bin/env python3

import sys
import os
import re
import argparse

parser = argparse.ArgumentParser(
    description='Analyse gtf output file for GALBA (or BRAKER) to generate a list of alternative transcripts required by OMARK (printed to STDOUT)')
parser.add_argument('-f', '--gtf_file', type=str, required=True,
                    help="Input gtf file (from GALBA or BRAKER)")
args = parser.parse_args()

# read gtf file into a dictionary that has gene name as key and a list of transcript names as value
gene_dict = {}
try:
    with open(args.gtf_file, "r") as file:
        for line in file:
            if re.search(r'gene_id \"([^"]+)\"', line):
                gid = re.search(r'gene_id \"([^"]+)\"', line).group(1)
            if re.search('transcript_id \"([^"]+)\"', line):
                tid = re.search('transcript_id \"([^"]+)\"', line).group(1)
            if re.search(r'transcript_id', line) and re.search(r'gene_id', line):
                if gid not in gene_dict:
                    gene_dict[gid] = {}
                if tid not in gene_dict[gid]:
                    gene_dict[gid][tid] = 1
except IOError:
    print("Could not read file:", args.gtf_file)
    sys.exit(1)

# print output to STDOUT in Omark specific format
for gid in gene_dict:
    first_tid = True
    tx_string = ""
    for tid in gene_dict[gid]:
        # tid.replace(" ", "")
        if first_tid:
            tx_string = tx_string + tid
            first_tid = False
        else:
            tx_string = tx_string + ";" + tid
    print(tx_string)
