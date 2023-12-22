#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(
    description='Check gtf file for errors and fix them if possible.')
parser.add_argument('-f', '--gtf_file', type=str, required=True,
                    help='Input gtf file (from GALBA because Pygustus ' +
                    'has a bug that we do not know how to fix, yet).')
parser.add_argument('-o', '--output', type=str, required=True,
                    help="Output gtf file.")
args = parser.parse_args()


def read_gtf(gtf_file):
    """
    Reads a GTF file and extracts CDS (Coding Sequence) information for each transcript.

    This function opens a GTF file and parses each line to find CDS features. For each CDS feature, 
    it extracts details like sequence ID, start and end positions, strand, frame, and transcript ID. 
    These details are organized and stored in dictionaries for each transcript. The function creates 
    three dictionaries: 'cds', and 'tx2str'. 'cds' is a nested dictionary where the first 
    key is the sequence ID and the second key is the transcript ID, containing a list of CDS feature 
    details. 'tx2str' maps 
    transcript IDs to strand information.

    Parameters:
    gtf_file (str): The file path of the GTF file to be read and processed.

    Returns:
    tuple: Returns a tuple containing three dictionaries:
        - cds (dict): A nested dictionary with CDS details for each transcript.
        - tx2str (dict): A dictionary mapping transcript IDs to strand information.

    Raises:
    IOError: If the GTF file cannot be opened or read, an IOError is raised with an error message.

    Note:
    The function expects the GTF file to be formatted correctly, particularly the CDS lines. 
    Incorrectly formatted files may lead to unexpected results or failures.
    """
    cds = {}
    try:
        with open(gtf_file, "r") as gtf_handle:
            for line in gtf_handle:
                if re.search(
                        r'\S+\t[^\t]+\tCDS\t\d+\t\d+\t\S+\t\S+\t\d\t.*transcript_id \"([^"]+)\"', line):
                    seq_id, st, en, stx, fr, tx_id = re.search(
                        r'(\S+)\t[^\t]+\tCDS\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\d)\t.*transcript_id \"([^"]+)\"', line).groups()
                    tx_id = re.sub(r'\"(\S+)\"', r'\1', tx_id)
                    tx_id = re.sub(r';', r'', tx_id)
                    if seq_id not in cds:
                        cds[seq_id] = {}
                    if tx_id not in cds[seq_id]:
                        cds[seq_id][tx_id] = []
                    cds[seq_id][tx_id].append(
                        {'start': int(st), 'end': int(en), 'strand': stx,
                         'frame': int(fr)})
    except IOError:
        print("Error: Failed to open file " + gtf_file + "!")
        exit(1)
    return cds
def check_gtf(cds_dict):
    """
    Checks the consistency and integrity of CDS (Coding Sequence) information within a GTF dataset.

    This function performs two main checks on the CDS data contained in 'cds_dict':
    1. Strand Consistency: Ensures that all CDS features for a given transcript (tx_id) 
       are on the same strand. If inconsistencies are found, the transcript ID is added to a list.
    2. Overlapping Features: Checks for overlapping CDS features within the same transcript. 
       If overlaps are found, the transcript ID is added to the list.

    The function returns a list of transcript IDs ('bad_lst') that failed these checks, indicating 
    potential issues in the GTF dataset.

    Parameters:
    cds_dict (dict): A nested dictionary containing CDS information. The first level of keys 
                     is sequence IDs (seq_id), and the second level of keys is transcript IDs 
                     (tx_id). Each entry contains a list of dictionaries with CDS feature details.

    Returns:
    list: A list of transcript IDs that failed the strand consistency or overlap checks.

    Note:
    The function assumes that the input 'cds_dict' is structured correctly, with appropriate 
    keys and values. Incorrectly structured input may lead to unexpected results.
    """
    bad_lst = []
    for seq_id in cds_dict:
        # check strand consistency
        for tx_id in cds_dict[seq_id]:
            strand = cds_dict[seq_id][tx_id][0]['strand']
            for cds in cds_dict[seq_id][tx_id]:
                if cds['strand'] != strand:
                    bad_lst.append(tx_id)
        # check whether two CDS features overlap
        for tx_id in cds_dict[seq_id]:
            for i in range(len(cds_dict[seq_id][tx_id])):
                for j in range(i+1, len(cds_dict[seq_id][tx_id])):
                    if cds_dict[seq_id][tx_id][i]['start'] < cds_dict[seq_id][tx_id][j]['end'] and cds_dict[seq_id][tx_id][i]['end'] > cds_dict[seq_id][tx_id][j]['start']:
                        bad_lst.append(tx_id)
    return bad_lst

def write_output(bad_lst, file):
    try:
        with open(file, "w") as file_handle:
            for tx_id in bad_lst:
                file_handle.write(tx_id + "\n")
    except IOError:
        print("Error: Failed to write output file " + file + "!")
        exit(1)


def main():
    # read gtf file
    gtf_dict = read_gtf(args.gtf_file)
    # check for errors
    bad_lst = check_gtf(gtf_dict)
    # write output list file
    write_output(bad_lst, args.output)

main()