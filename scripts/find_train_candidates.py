#!/usr/bin/env python3

import argparse
import re
from Bio import SeqIO

__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2023. All rights reserved."
__license__ = "Artistic Licsense"
__credits__ = "ChatGPT"
__version__ = "1.0.0"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "development"


parser = argparse.ArgumentParser(
    description='Parse miniprot gtf output, select the highest scoring transcript per gene range on the same strand.')
parser.add_argument('-m', '--miniprot_gtf', required=True, type=str,
                    help='Miniprot output in gtf format.')
parser.add_argument('-o', '--output', required=True, type=str,
                    help='Output file with training gene candidates in gtf format.')
args = parser.parse_args()


def read_gtf(filename):
    """
    Reads a GTF file and returns a dictionary where each key is a chromosome and each value
    is a list of genes on that chromosome, where each gene is represented as a dictionary with
    the gene's attributes as keys and their values as values.
    """
    chromosomes = {}
    with open(filename, "r") as gtf_file:
        for line in gtf_file:
            if line.startswith("#") or re.search(r'\*', line):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "gene":
                continue
            attributes = dict(item.strip().split(" ") for item in fields[8].split(";") if item.strip())
            chromosome = fields[0]
            strand = fields[6]
            start = int(fields[3])
            end = int(fields[4])
            score = int(fields[5])
            gene_id = attributes["gene_id"]
            if chromosome not in chromosomes:
                chromosomes[chromosome] = []
            chromosomes[chromosome].append({
                "gene_id": gene_id,
                "start": start,
                "end": end,
                "strand": strand,
                "score": score
            })
    return chromosomes


def identify_highest_score_genes(chromosomes):
    """
    Identifies the gene with the highest score for each cluster of overlapping genes on the same strand.
    Returns a dictionary where each key is a chromosome and each value is a list of genes, where each
    gene is represented as a dictionary with the gene's attributes as keys and their values as values.
    """
    highest_score_genes = {}
    for chromosome, genes in chromosomes.items():
        sorted_genes = sorted(genes, key=lambda x: x["start"])
        current_cluster = []
        for gene in sorted_genes:
            if not current_cluster:
                current_cluster.append(gene)
            elif gene["start"] <= current_cluster[-1]["end"]:
                current_cluster.append(gene)
            else:
                highest_score_gene = max(current_cluster, key=lambda x: x["score"])
                if chromosome not in highest_score_genes:
                    highest_score_genes[chromosome] = []
                highest_score_genes[chromosome].append(highest_score_gene)
                current_cluster = [gene]
        if current_cluster:
            highest_score_gene = max(current_cluster, key=lambda x: x["score"])
            if chromosome not in highest_score_genes:
                highest_score_genes[chromosome] = []
            highest_score_genes[chromosome].append(highest_score_gene)
    return highest_score_genes


def convert_to_gene_id_dict(highest_score_genes):
    """
    Takes the output of identify_highest_score_genes() and converts it into a dictionary
    where the gene ID is the key and the value is an empty string.
    """
    gene_id_dict = {}
    for chromosome, genes in highest_score_genes.items():
        for gene in genes:
            gene_id = gene["gene_id"].replace('"', '')
            gene_id_dict[gene_id] = ""
    return gene_id_dict


def print_genes_from_gtf(gtf_file, gene_id_dict, out_file):
    """
    Reads a GTF file and prints all lines corresponding to genes whose IDs are stored in the
    given dictionary.
    """
    try:
        with open(gtf_file) as f:
            try:
                with open(out_file, "w") as o:
                    for line in f:
                        if line.startswith("#") or re.search(r'\*', line):
                            # skip comments and output lines with stars
                            continue
                        fields = line.split("\t")
                        attributes = dict(item.strip().split(" ") for item in fields[8].split(";") if item.strip())
                        gene_id = attributes["gene_id"].strip('"')
                        if gene_id in gene_id_dict:
                            # print lines for matching gene IDs
                            o.write(line)
            except IOError:
                print("Error: failed to open file " + out_file + " for writing!")
    except IOError:
        print("Error: failed to open file " + gtf_file + "for reading!")
        exit(1)


# read gtf file
all_genes = read_gtf(args.miniprot_gtf)

# find the genes with largest score per cluster over overlapping genes
best_genes = identify_highest_score_genes(all_genes)

# make a compact dict for easy lookup when reading gtf
selected_gids = convert_to_gene_id_dict(best_genes)

# print the selected genes with maximal score
print_genes_from_gtf(args.miniprot_gtf, selected_gids, args.output)
