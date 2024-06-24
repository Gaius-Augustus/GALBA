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
    
        result = subprocess.run(hisat2_build_command, capture_output=True)
    
        if result.returncode == 0:
            print("Indexing completed successfully.")
        else:
            print("Error during indexing:")
            print(result.stderr)

    except Exception as e:
        print("An error occurred: {e}") 

def mapping(output, reads_file):
    try:
        hisat2_command = [
            "hisat2",
            "-f",                        
            "-x", output,            
            "-U", reads_file,            
            "-S", "output.sam",            
            "--no-spliced-alignment"        
        ]
        
        print("Running command:", " ".join(hisat2_command))
        
        result = subprocess.run(hisat2_command)
        
        if result.returncode == 0:
            print("Mapping completed successfully")
        else:
            print("Error during mapping:")
            print(result.stderr)
    
    except Exception as e:
        print("ERROR") 

def sam_to_bam(samFile):
    try:
        samtools_command = [
            "samtools",
            "view",
            "-b",
            samFile,
            "-o",
            "output.bam"
        ]

        print("Running command:", "".join(samtools_command))

        result = subprocess.run(samtools_command)

        if result.returncode == 0:
            print("Conversion from .sam to .bam file completed successfully")

        else:
            print("Error during conversion:")
            print(result.stderr)

    except Exception:
        print("ERROR samtools")

#parser = argparse.ArgumentParser(description='Genome annotation pipeline')

#Input: Genome file and reads file and defining output name
genome_fasta = sys.argv[1]
output = sys.argv[2]
reads_file = sys.argv[3]

if not os.path.isfile(genome_fasta):
    print("Error: The file {genome_fasta} does not exist.")

else:
    indexing(genome_fasta, output)

if not os.path.isfile(reads_file):
    print("Error: The file {reads_file} does not exist.")

else:
    mapping(output, reads_file)
    sam_to_bam("output.sam")
        