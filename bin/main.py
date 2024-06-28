#!/usr/bin/env python3

#import argparse
import subprocess
import sys
import os

def indexing(genome_fasta, output):
    try:
        hisat2_build_command = [
            "hisat2-build",         
            genome_fasta,           
            output         
        ]
    
        print("Running command:", " ".join(hisat2_build_command))
    
        result = subprocess.run(hisat2_build_command, capture_output=False)
    
        if result.returncode == 0:
            print("Indexing completed successfully.")
        else:
            print("Error during indexing:")
            print(result.stderr)

    except Exception as e:
        print("An error occurred: {e}") 

def mapping(output_indexing, reads_file, output_sam):
    try:
        hisat2_command = [
            "hisat2",
            "-f",                        
            "-x", output_indexing,            
            "-U", reads_file,            
            "-S", output_sam,            
            "--no-spliced-alignment"        
        ]
        
        print("Running command:", " ".join(hisat2_command))
        
        result = subprocess.run(hisat2_command, capture_output=True)
        
        if result.returncode == 0:
            print("Mapping completed successfully")
        else:
            print("Error during mapping:")
            print(result.stderr)
    
    except Exception as e:
        print("ERROR") 

def sam_to_bam(samFile, output_bam):
    try:
        samtools_command = [
            "samtools",
            "sort",
            samFile,
            "-o",
            output_bam
        ]

        print("Running command:", "".join(samtools_command))

        result = subprocess.run(samtools_command, capture_output=True)

        if result.returncode == 0:
            print("Conversion from .sam to .bam file completed successfully")

        else:
            print("Error during conversion:")
            print(result.stderr)

    except Exception:
        print("ERROR samtools")

def fastq_to_fasta(fastqFile, output_fasta):
    stringtie_path = os.path.expanduser("/home/s-amknut/GALBA/tools/seqtk/./seqtk")

    seqtk_command = [
        stringtie_path,
        "seq",
        "-a",
        fastqFile,
        ">",
        output_fasta
    ]

    print("Running command:", " ".join(seqtk_command))

    result = subprocess.run(seqtk_command, capture_output=True)

    if result.returncode == 0:
        print("Conversion from .fastq to .fasta file completed successfully")

    else:
        print("Error in converting .fastq to .fasta file: ")
        print(result.stderr)
#parser = argparse.ArgumentParser(description='Genome annotation pipeline')

#Input: Genome file and reads file and defining output name
genome_fasta = sys.argv[1]
output_indexing = sys.argv[2]
reads_file = sys.argv[3]

if not os.path.isfile(genome_fasta):
    print("Error: The file {genome_fasta} does not exist.")

else:
    indexing(genome_fasta, output_indexing)

if not os.path.isfile(reads_file):
    print("Error: The file {reads_file} does not exist.")

else:
    reads_fasta = "reads.fasta"
    fastq_to_fasta(reads_file, reads_fasta)
    output_sam = "output.sam"
    mapping(output_indexing, reads_fasta, output_sam)
    sam_to_bam(output_sam, "output.bam") 
        