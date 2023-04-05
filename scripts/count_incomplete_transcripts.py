#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(
    description='Analyse gff output of AUGUSTUS with respect to complete and incomplete predictions')
parser.add_argument('-f', '--gff_file', type=str, required=True,
                    help="Input gff file (from AUGUSTUS)")
args = parser.parse_args()

# Define a function to extract the transcript ID from a GTF line                                                                                              
def extract_transcript_id(line):
    if re.search(r'transcript_id \"([^"]+)\"', line):
        return re.search(r'transcript_id \"([^"]+)\"', line).group(1)
    else:
        return None

# Define sets to keep track of all transcripts and transcripts with both start and stop codons                                                                
transcripts = {}

# First, read in the file
with open(args.gff_file, "r") as f:
    for line in f:
        if re.search(r'\tAUGUSTUS\t', line):
            fields = line.strip().split("\t")
            feature = fields[2]
            if feature == "CDS":
                transcript_id = extract_transcript_id(line)
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = 0
            elif feature == "start_codon":
                transcript_id = extract_transcript_id(line)
                if transcript_id in transcripts:
                    transcripts[transcript_id] = transcripts[transcript_id] + 1
                else:
                    transcripts[transcript_id] = 1
            elif feature == "stop_codon":
                transcript_id = extract_transcript_id(line)
                if transcript_id in transcripts:
                    transcripts[transcript_id] = transcripts[transcript_id] + 1
                else:
                    transcripts[transcript_id] = 1

# Print the results
print("Total transcripts:", len(transcripts.keys()))
completes = 0
for tx in transcripts.keys():
#    print(tx)
#    print(transcripts[tx])
    if transcripts[tx] >= 2:
        completes = completes + 1

print("Transcripts with both start and stop codons:", completes)

