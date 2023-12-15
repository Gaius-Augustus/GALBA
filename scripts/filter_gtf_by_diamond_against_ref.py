#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import re

__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2023. All rights reserved."
__credits__ = ""
__license__ = "Artistic License 2.0"
__version__ = "1.0.0"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "development"

argparser = argparse.ArgumentParser(description = 'Filter genes predicted by ' +
                                    'AUGUSTUS within e.g. GALBA by simply ' +
                                    'discarding all genes that are not ' +
                                    'supported by a DIAMOND hit to the ' +
                                    'reference protein set. Warning: this is ' +
                                    'a poor approach if you do not have a ' +
                                    'high quality reference protein set from ' +
                                    'close relatives. This approach may lead ' +
                                    'to discarding species-specific gene ' +
                                    'models.')
argparser.add_argument('-r', '--ref_prots', type=str, required = True,
                            help = 'FASTA file with reference protien sequences')
argparser.add_argument('-g', '--gtf', type=str, required = False,
                            help = 'gtf file with gene models to be filtered')
argparser.add_argument('-a', '--pred_prots', type=str, required = False,
                            help = 'FASTA file with predicted protein sequences')
argparser.add_argument('-o', '--out_gtf', type=str, required = True,
                            help = 'Output file in gtf format')
argparser.add_argument('-t', '--threads', type=str, required = False, default = 1,
                       help = 'Number of threads to use for running DIAMOND')
argparser.add_argument('-d', '--diamond_dir', type=str, required = True,
                            help = 'Directory that contains DIAMOND binaries')

args = argparser.parse_args()

def run_diamond(ref_prots, pred_prots, threads, diamond_dir):
    """
    Run DIAMOND to search for homologous proteins in the reference protein set
    """
    diamond_out_file = "diamond_out.tsv"
    diamond_out = open(diamond_out_file, "w")
    diamond_cmd = diamond_dir + "/diamond blastp -d " + ref_prots + " -q " + pred_prots + " -o " + diamond_out_file + " -p " + threads + " --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    diamond_out.write("Running DIAMOND with command: " + diamond_cmd + "\n")
    diamond_out.flush() 
    diamond_out.close()
    subprocess.call(diamond_cmd, shell=True)
    return diamond_out_file

def filter_gtf(gtf, diamond_out_file, out_gtf):
    """
    Read DIAMOND output, identify transcript names that have at least one hit,
    read gtf file, discard all genes that do not have a hit,
    write filtered gtf file
    """
    # reading DIAMOND output
    try:
        with open(diamond_out_file, "r") as diamond_out:
            diamond_hits_tx = {}
            diamond_hits_gene = {}
            for line in diamond_out:
                line = line.rstrip()
                fields = line.split("\t")
                query = fields[0]
                if query not in diamond_hits_tx:
                    diamond_hits_tx[query] = 1
                    match = re.search(r'(g\d+)\.t\d+', query)
                    diamond_hits_gene[match.group(1)] = 1
    except IOError:
        print("ERROR: DIAMOND output file " + diamond_out_file + " could not be opened")
        sys.exit(1)

    # reading gtf file to filter out genes that do not have a hit
    try:
        with open(gtf, "r") as gtf_in:
            try:
                with open(out_gtf, "w") as gtf_out:
                    for line in gtf_in:
                        line = line.rstrip()
                        fields = line.split("\t")
                        if fields[2] == "gene":
                            if fields[8] in diamond_hits_gene:
                                gtf_out.write(line + "\n")
                        elif fields[2] == "transcript":
                            if fields[8] in diamond_hits_tx:
                                gtf_out.write(line + "\n")
                        else:
                            match = re.search('transcript_id "([^"]+)";', fields[8])
                            txid = match.groups()[0]
                            if txid in diamond_hits_tx:
                                gtf_out.write(line + "\n")
                            else:
                                print("Discarding " + txid)
            except IOError:
                print("ERROR: gtf file " + out_gtf + " could not be opened")
                sys.exit(1)
    except IOError:
        print("ERROR: gtf file " + gtf + " could not be opened")
        sys.exit(1)
        

def main():
    diamond_out_file = run_diamond(args.ref_prots, args.pred_prots, args.threads, args.diamond_dir)
    filter_gtf(args.gtf, diamond_out_file, args.out_gtf)
    os.remove(diamond_out_file)

main()