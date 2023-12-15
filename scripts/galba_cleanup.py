#!/usr/bin/env python3

# ==============================================================
# Katharina Hoff
#
# ==============================================================

import os
import shutil
import argparse

# Define the list of files to delete
files_to_delete = [
    "firsttest.stdout", "genome.fa", "getAnnoFasta.augustus.ab_initio.stdout",
    "getAnnoFasta.augustus.hints.stdout", "secondtest.stdout",
    "train.gb",
    "aug_hints.lst", "aa2nonred.stdout", "augustus.hints.tmp.gtf", "fourthtest.stdout", "gbFilterEtraining.stdout",
    "genes.gtf", "genes_in_gb.gtf",
    "singlecds.hints", "stops.and.starts.gff", "train.gb.test", "train.gb.train", "train.gb.train.test",
    "train.gb.train.train", "traingenes.good.fa", "downsample_traingenes.log", "firstetraining.stdout",
    "secondetraining.stdout", "startAlign_gth.log", "protein_alignment_gth.gff3", "ex1.cfg", "getAnnoFastaFromJoingenes.augustus.hints.stdout",
    "genome.fa.cidx", "getAnnoFastaFromJoingenes.augustus.hints_tmp.stdout", "getAnnoFastaFromJoingenes.augustus.ab_initio_tmp.stdout",
    "augustus.ab_initio.tmp.gtf", "augustus.ab_initio.gff", "augustus.hints.tmp.gtf",
    "getAnnoFastaFromJoingenes.augustus.hints_hints.stdout", "getAnnoFastaFromJoingenes.augustus.ab_initio_.stdout",
    "getAnnoFastaFromJoingenes.augustus.hints_.stdout", "startAlign.stdout", "cmd.log", "etrain.bad.lst", "gene_stat.yaml",
    "good_genes.lst", "nonred.loci.lst", "nuc.fasta", "proteins.fa",
    "train.f.gb", "traingenes.good.gtf", "traingenes.good.nr.fa", "uniqueSeeds.gtf", "genome.mpi", "pygustus_hints.out", "pygustus_hints.py",
    "gff2gbSmallDNA.stderr", "miniprot_representatives.gff", "protein_alignment_miniprot.aln", "protein_alignment_miniprot.gff", "hc.gff",
    "miniprothint.gff", "thirdtest.stdout", "miniprot_trainingGenes.gff"
]

def main():
    parser = argparse.ArgumentParser(description="Delete galba.pl output files that are usually not required for downstream analysis")
    parser.add_argument("--wdir", required=True, help="Output directory of galba.pl job")
    args = parser.parse_args()

    # Check if the specified directory exists
    if not os.path.exists(args.wdir):
        print(f"ERROR: Directory {args.wdir} does not exist!")
        return

    # Delete individual files
    for filename in files_to_delete:
        file_path = os.path.join(args.wdir, filename)
        if os.path.exists(file_path):
            print(f"Deleting file {file_path}")
            os.remove(file_path)

    # Delete directories with names like "align_gth0", "align_gth1", and so on
    gth_index = 0
    while os.path.exists(os.path.join(args.wdir, f"align_gth{gth_index}")):
        dir_path = os.path.join(args.wdir, f"align_gth{gth_index}")
        print(f"Deleting directory {dir_path}")
        shutil.rmtree(dir_path)
        gth_index += 1

    # Create a new directory for archived files
    archive_dir = os.path.join(args.wdir, "archive")
    if not os.path.exists(archive_dir):
        os.mkdir(archive_dir)

    # Move specific files to the archive directory
    archive_files = [
        "augustus.hints.gff", "train2.gb.train.test", "train2.gb.train.train",
        "train2.gb.test", "train2.gb.train", "train2.gb", "miniprot_representatives.gtf", "miniprot.gff",
        "filter_gtf_by_diamond_against_ref.stdout"
    ]
    for filename in archive_files:
        file_path = os.path.join(args.wdir, filename)
        if os.path.exists(file_path):
            print(f"Moving file {file_path} to archive directory")
            shutil.move(file_path, os.path.join(archive_dir, filename))

if __name__ == "__main__":
    main()
