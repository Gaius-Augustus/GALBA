#!/usr/bin/env python3
# ==============================================================
# Katharina Hoff
#
# ==============================================================

import sys
import os

def find_strand(seqname, start, end, allowed):
    reverse = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    seq = annos.get(seqname, '')
    if not seq:
        print(f"WARNING: '{seqname}' does not match any sequence in the fasta file. Maybe the two files do not belong together.", file=sys.stderr)
        return '.'

    intron_type = seq[start - 1:start + 1].lower() + seq[end - 2:end].lower()
    reverse_intron_type = ''.join(reverse[base] for base in reversed(intron_type))

    has_plus = any(intron_type == allowed_type for allowed_type in allowed)
    has_minus = any(reverse_intron_type == allowed_type for allowed_type in allowed)

    if has_plus and has_minus:
        return 'b'
    elif has_minus:
        return '-'
    elif has_plus:
        return '+'
    else:
        return '.'

def get_score(column):
    mult_match = re.search(r'mult=(\d+)', column)
    if mult_match:
        return mult_match.group(1)
    else:
        return '1'

if len(sys.argv) < 3:
    print("Usage: filterIntronsFindStrand.py genome.fa introns.gff [OPTIONS] > introns.s.f.gff")
    sys.exit(1)

genome = sys.argv[1]
introns = sys.argv[2]

# Set allowed splice site types
allowed = ["gtag", "gcag", "atac"]

mult_score = False
annos = {}
seqname = None
seq = ''

with open(genome, 'r') as fasta_file:
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            if seqname is not None:
                annos[seqname] = seq
            seqname = line[1:]
            seq = ''
        else:
            seq += line
    annos[seqname] = seq

with open(introns, 'r') as introns_file:
    for line in introns_file:
        line = line.strip().split('\t')
        if len(line) == 9:
            if '.' in line[6]:
                strand = find_strand(line[0], int(line[3]), int(line[4]), allowed)
                if strand in ['+', '-', 'b']:
                    score = get_score(line[8]) if mult_score else line[5]
                    if strand == 'b':
                        print(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}\t{score}\t+\t{line[7]}\t{line[8]}")
                        print(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}\t{score}\t-\t{line[7]}\t{line[8]}")
                    else:
                        print(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}\t{score}\t{strand}\t{line[7]}\t{line[8]}")
            else:
                print('\t'.join(line))

