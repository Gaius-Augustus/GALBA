#!/usr/bin/env python

import argparse
import subprocess
import sys
import os

def fastq_to_fasta(fastqFile):
    seqtk_path = os.path.expanduser("/home/s-amknut/GALBA/tools/seqtk/./seqtk")
    output_fasta = os.path.expanduser("/home/s-amknut/GALBA/bin/reads.fa")

    seqtk_command = [
        seqtk_path,
        "seq",
        "-a",
        fastqFile,
        ">",
        output_fasta
    ]

    print("Running command:", " ".join(seqtk_command))

    with open(output_fasta, 'w') as output_file:

        result = subprocess.run(seqtk_command, stdout=output_file, capture_output=False)

    if result.returncode == 0:
        print("Conversion from .fastq to .fasta file completed successfully")

    else:
        print("Error in converting .fastq to .fasta file: ")
        print(result.stderr)

def file_format(file):
    with open(file, 'r') as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()
        third_line = f.readline().strip()

        if first_line.startswith('>'):
            format = "fasta"
            print("Input file is fasta format")
            return format
        elif first_line.startswith('@') and third_line.startswith('+'):
            format = "fastq"
            print("Input file is fastq format")
            return format
        else:
            format = "unknown"
            print("File format is unknown")
            return format

def indexing(genome_fasta, output):
    if not os.path.isfile(genome_fasta):
        print("Error: The file {genome_fasta} does not exist.")

    else:
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
    if not os.path.isfile(reads_file):
        print("Error: The file {reads_file} does not exist.")

    else: 
        if file_format(reads_file) == "fastq":
            print("Converting fastq to fasta file")
            fastq_to_fasta(reads_file)
            reads_file = "reads.fa"
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

def assembling(bamFile):
    stringtie_path = os.path.expanduser("/home/s-amknut/GALBA/tools/stringtie/./stringtie")
    output_gtf = os.path.expanduser("/home/s-amknut/GALBA/bin/assembly.gtf")

    try:
        stringtie_command = [
            stringtie_path,
            bamFile,
            "-o",
            output_gtf
        ]
    
        print("Running command:", "".join(stringtie_command))

        result = subprocess.run(stringtie_command, capture_output=True)

        if result.returncode == 0:
            print("Assembled reads successfully")

        else:
            print("Error during Assembly:")
            print(result.stderr)

    except Exception:
        print("ERROR stringtie")

parser = argparse.ArgumentParser()  
parser.add_argument('--g', help='Genome file', required=True)
parser.add_argument('--r', help='Reads file', required=True)
#parser.add_argument('--o', help='Output name', required=True)

args = parser.parse_args()
genome_file = args.g
reads_file = args.r
output_indexing = "output_indexing"

indexing(genome_file, output_indexing)
mapping(output_indexing, reads_file, "output.sam")
sam_to_bam("output.sam", "output.bam") 
assembling("output.bam")


"""
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
"""      