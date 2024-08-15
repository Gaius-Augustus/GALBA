#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(
    description='Analyse gtf output file for GALBA (or BRAKER) with respect of CDS exons per transcript')
parser.add_argument('-f', '--gtf_file', type=str, required=True,
                    help="Input gtf file (from GALBA or BRAKER)")
parser.add_argument('-o', '--intron_out', type=str, required=False,
                    help="output text file with number of introns per transcripts")
args = parser.parse_args()

# read gtf file into a dictionary that has transcript name as key and number of exons as value
tx_dict = {}
try:
    with open(args.gtf_file, "r") as gtf_file:
        for line in gtf_file:
            if re.search(r'\tCDS\t', line):
                if re.search(r'transcript_id \"([^"]+)\"', line).group(1) not in tx_dict:
                    if re.search(r'\tAUGUSTUS\t', line):
                        # ignore alternative transcripts because they do not compare to single exon transcripts, take arbitrarily the first one
                        if re.search(r't1', re.search(r'transcript_id \"([^"]+)\"', line).group(1)):
                            tx_dict[re.search(r'transcript_id \"([^"]+)\"', line).group(1)] = 1
                    else:
                        # GeneMark should not have alternatives
                        if re.search(r'transcript_id \"([^"]+)\"', line).group(1):
                            tx_dict[re.search(r'transcript_id \"([^"]+)\"', line).group(1)] = 1
                else:
                    if re.search(r'\tAUGUSTUS\t', line):
                        if re.search(r't1', re.search(r'transcript_id \"([^"]+)\"', line).group(1)):
                            tx_dict[re.search(r'transcript_id \"([^"]+)\"', line).group(1)] += 1
                    else:
                        # GeneMark should not have alternatives
                        if re.search(r'transcript_id \"([^"]+)\"', line).group(1):
                            tx_dict[re.search(r'transcript_id \"([^"]+)\"', line).group(1)] += 1
except IOError:
    print("Could not read file:", args.gtf_file)
    sys.exit(1)

# output some information about the dictionary
print("Number of transcripts:", len(tx_dict))
print("Largest number of exons in all transcripts:", max(tx_dict.values()))
print("Monoexonic transcripts:", list(tx_dict.values()).count(1))
print("Multiexonic transcripts:", len(tx_dict) - list(tx_dict.values()).count(1))
print("Mono:Mult Ratio:", round(list(tx_dict.values()).count(1) / (len(tx_dict) - list(tx_dict.values()).count(1)), 2))

# print an ASCII boxplot of data in tx_dict
# get the largest number of exons in all transcripts
max_exons = max(tx_dict.values())
# get the smallest number of exons in all transcripts
min_exons = min(tx_dict.values())
# get the 25% quartile of values in tx_dict
q25 = sorted(tx_dict.values())[int(len(tx_dict) * 0.25)]
# get the 50% quartile of values in tx_dict
q50 = sorted(tx_dict.values())[int(len(tx_dict) * 0.5)]
# get the 75% quartile of values in tx_dict
q75 = sorted(tx_dict.values())[int(len(tx_dict) * 0.75)]


# print the boxplot
print("Boxplot of number of exons per transcript:")
print("Min:", min_exons)
print("25%:", q25)
print("50%:", q50)
print("75%:", q75)
print("Max:", max_exons)

if args.intron_out:
    with open(args.intron_out, 'w') as out_file:
        for value in tx_dict.values():
            out_file.write(str(int(value)-1) + "\n")
        
