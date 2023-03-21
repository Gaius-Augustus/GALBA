#!/usr/bin/env python3

# take a multiple fasta amino acid sequence file
# and output a file that is non-redundant:
#
# the percent identity value between each pair of different sequence is below a threshold
# and each sequence name occurrs only once
# we do not select randomly, but keep the longest sequence

import argparse
import re
import os
import subprocess

__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2023. All rights reserved."
__license__ = "Artistic Licsense"
__credits__ = "ChatGPT"
__version__ = "1.0.0"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "development"

parser = argparse.ArgumentParser(
    description='make a protein file non-redundant.')
parser.add_argument('-i', '--input_file', required=True, type=str,
                    help='Protein sequence input file in multiple fasta format.')
parser.add_argument('-o', '--output', required=True, type=str,
                    help='Output file with non-redundant protein sequences.')
parser.add_argument('-m', '--maxid', required=False, type=float, default=0.8,
                    help='Maximum percent identity between to sequences (#identical aa) / (length of shorter sequence) default: 0.8')
parser.add_argument('-d', '--diamond_path', required=False, type=str, default="diamond",
                    help='Path to diamond')
parser.add_argument('-t', '--threads', required=False, type=int, default=1,
                    help='Number of threads to be used by DIAMOND')
args = parser.parse_args()

'''set_DIAMOND_PATH() sets the path to diamond executable'''
def set_DIAMOND_PATH():
    if 'DIAMOND_PATH' in os.environ:
        if os.path.exists(os.environ['DIAMOND_PATH']):
            return os.environ['DIAMOND_PATH']
    if args.diamond_path:
        if os.path.isdir(args.diamond_path):
            return diamond_path
    try:
        epath = subprocess.check_output(['which', 'diamond']).decode().strip()
        if os.path.isdir(os.path.dirname(epath)):
            return os.path.dirname(epath)
    except subprocess.CalledProcessError:
        print("Error: DIAMOND not found. Please install DIAMOND and add it to your PATH or set the environment variable DIAMOND_PATH to the DIAMOND executable.")
        sys.exit(1)

diamond_bin = set_DIAMOND_PATH()
