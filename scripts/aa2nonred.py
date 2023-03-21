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
import sys
from inspect import currentframe, getframeinfo

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

''' create_seqname_lst(in.fa, out.lst) creates a list of unique sequence names from a fasta file '''

def create_seqname_lst(infa, outlst):
    seqnames = {}
    try:
        with open(infa, "r") as f:
            with open(outlst, "w") as g:
                for line in f:
                    if line.startswith(">"):
                        if re.search(r'>(\S+)', line).group(1) not in seqnames:
                            seqnames[re.search(r'>(\S+)', line).group(1)] = 1
    except IOError:
        print("Error: cannot open file " + infa)
        sys.exit(1)
    try:
        with open(outlst, "w") as g:
            for seqname in seqnames:
                g.write(seqname + "\n")
    except IOError:
        print("Error: cannot open file " + outlst)
        sys.exit(1)

def run_simple_process(args_lst):
    try:
        print("Trying to execute the following command:")
        print(" ".join(args_lst))
        result = subprocess.run(
            args_lst, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Suceeded in executing command.")
        if(result.returncode == 0):
            return(result)
        else:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Return code of subprocess was " +
                  str(result.returncode) + str(result.args))
            quit(1)
    except subprocess.CalledProcessError as grepexc:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed executing: ",
              " ".join(grepexec.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        quit(1)

''' create DIAMOND database from fasta file and run reciprocal search '''

def run_diamond(diamond_bin, infa, outdaa):
    run_simple_process([diamond_bin, '--in', infa, '-d', infa + "diamond.db", '--threads', str(args.threads)])
    run_simple_process([diamond_bin, '--db', infa + "diamond.db", '--outfmt', '0', '--query', infa, '--threads', str(args.threads), '--out', outdaa])

diamond_bin = set_DIAMOND_PATH() + "/diamond"
create_seqname_lst(args.input_file, args.input_file + ".nonred.db")
run_diamond(diamond_bin, args.input_file + ".nonred.db", args.input_file + ".nonred.daa")
