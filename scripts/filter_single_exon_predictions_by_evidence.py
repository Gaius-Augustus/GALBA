#!/usr/bin/env python3

import argparse
import re
import sys

parser = argparse.ArgumentParser(
    description='Analyse gff output file AUGUSTUS for identify single exon transcripts that have zero support by evidence')
parser.add_argument('-f', '--gff_file', type=str, required=True,
                    help="Output of AUGUSTUS in gff format")
args = parser.parse_args()

# loop over gff_file
# search for lines with \tsingle\t and extract transcript ID
# after extraction, search for %of transcript supported by hints (any source): 0
# if found, print transcript ID
look = False
try:
    with open(args.gff_file, "r") as gff_file:
        for line in gff_file:
            if re.search(r'\tsingle\t', line):
                txid = re.search(r'transcript_id \"([^"]+)\"', line).group(1)
                look = True
            if look is True:
                if re.search(r'%of transcript supported by hints \(any source\): (\d+)', line):
                    support = re.search(r'%of transcript supported by hints \(any source\): (\d+)', line).group(1)
                    if support == '0':
                        print(txid)
                    look = False
except IOError:
    print("Could not read file:", args.gff_file)
    sys.exit(1)