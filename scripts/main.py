#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os

#VARIABLES 
cores = 4


#FUNCTIONS
def check_input(genome_file, reads_file):
    if not os.path.isfile(genome_file):
        print("Error: Could not find the genome file.")
        sys.exit(1)
    if not os.path.isfile(reads_file):
        print("Error: Could not find the reads file.")
        sys.exit(1)
    if file_format(reads_file) == "fastq":
        print("Input reads file is in fastq format")
        fastq_to_fasta(reads_file)
    if file_format(genome_file) == "fastq":
        print("Input genome file is in fastq format")
        fastq_to_fasta(genome_file)
    if file_format(reads_file) == "unknown":
        sys.exit("Error Reads file: Unknown file format")
    if file_format(genome_file) == "unknown":
        sys.exit("Error Genome file: Unknown file format")

def fastq_to_fasta(fastqFile):
    seqtk_path = os.path.expanduser("/home/s-amknut/GALBA/tools/seqtk/./seqtk")
    try:
        seqtk_command = [
            seqtk_path,
            "seq",
            "-a",
            fastqFile,
            "> reads.fa",
        ]

        print("Converting .fastq to .fasta file...")

        with open("reads.fa", 'w') as output_file:
            result = subprocess.run(seqtk_command, stdout=output_file, capture_output=False)

        if result.returncode == 0:
            print("Conversion from .fastq to .fasta file completed successfully")

        else:
            print("Error in converting .fastq to .fasta file: ")
            print(result.stderr)
            #sys.exit(1)
    
    except Exception:
        print("Could not run seqtk command.")
       # sys.exit(1)

def file_format(file):
    with open(file, 'r') as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()
        third_line = f.readline().strip()

        if first_line.startswith('>'):
            format = "fasta"
            return format
        elif first_line.startswith('@') and third_line.startswith('+'):
            format = "fastq"
            return format
        else:
            format = "unknown"
            return format

def indexing(genome_fasta, output):
    hisat2_build_path = os.path.expanduser("/home/s-amknut/GALBA/tools/hisat2/hisat2-build")
    try:
        hisat2_build_command = [
            hisat2_build_path, 
            "--quiet",
            "-p",
            str(cores),        
            genome_fasta,           
            output         
        ]

        print("Building genome index...")
    
        result = subprocess.run(hisat2_build_command, capture_output=True)
    
        if result.returncode == 0:
            print("Indexing completed successfully.")
        else:
            print("Error during indexing:")
            print(result.stderr)

    except Exception:
        print("Could not run hisat2-build command.") 

def mapping(output_indexing, reads_file, output_sam):
    hisat2_path = os.path.expanduser("/opt/hisat2/hisat2")
    if file_format(reads_file) == "fastq":
        fastq_to_fasta(reads_file)
        reads_file = "reads.fa"
    try:
        #which Pfad finden
        hisat2_command = [
            "hisat2",
            "-f",                        
            "-x", output_indexing,            
            "-U", reads_file,            
            "-S", output_sam,            
            "--no-spliced-alignment"        
        ]
        
        print("Mapping reads to genome...")
        
        result = subprocess.run(hisat2_command, capture_output=True)
        
        if result.returncode == 0:
            print("Mapping completed successfully")
        else:
            print("Error during mapping:")
            print(result.stderr)
    
    except Exception as e:
        print("ERROR") 

def sam_to_bam(samFile, output_bam):
    if not os.path.isfile(samFile):
        print("Error: The bamfile does not exist.")

    try:
        samtools_command = [
            "samtools",
            "sort",
            samFile,
            "-o",
            output_bam
        ]

        #print("Running command:", "".join(samtools_command))
        print("Converting .sam to .bam file...")
        result = subprocess.run(samtools_command, capture_output=True)

        if result.returncode == 0:
            print("Conversion from .sam to .bam file completed successfully")

        else:
            print("Error during conversion:")
            print(result.stderr)

    except Exception:
        print("Could not run samtools command.")

def assembling(bamFile, output_gtf):
    if not os.path.isfile(bamFile):
        print("Error: The bamfile does not exist.")

    try:
        print("Assembling the reads...")
        command = "/home/s-amknut/GALBA/tools/stringtie2/stringtie -p "+str(cores)+" -o "+output_gtf+" "+bamFile
        result = os.system(command)
        if result== 0:
            print("Assembled reads successfully")

        else:
            print("Error during Assembly")

    except Exception:
        print("Could not run stringtie command.")

def orfsearching(assembly_gtf, genome_fa):
    try:
        print("Preparing Transcripts using gffread...")
        gffread_command = "/home/s-amknut/GALBA/tools/gffread/./gffread " + assembly_gtf+" "+ genome_fa +" > transcripts.fasta"
        print("Running command:", "".join(gffread_command))
        result = os.system(gffread_command)
        #result = subprocess.run(gffread_command, capture_output=True)
        if result.returncode == 0:
            print("Transcripts prepared successfully")
        else:
            print("Error during preparing transcripts")
            print(result.stderr)    
    
    except Exception:
        print("Could not run gffread command.")

    try:
        print("Searching for ORFs...")
        transdecoder_command = [
            "/home/s-amknut/GALBA/tools/eviann/.TransDecoder.LongOrfs",
            " -t transcripts.fasta"
        ]

        result= subprocess.run(transdecoder_command, capture_output=True)
        if result.returncode == 0:
            print("ORFs found successfully")
        else:
            print("Error during ORF search")
            print(result.stderr)
    except Exception:
        print("Could not run TransDecoder command.")


#MAIN
parser = argparse.ArgumentParser()  
parser.add_argument('-g', help='Genome file', required=True)
parser.add_argument('-r', help='Reads file', required=True)

args = parser.parse_args()
genome_file = args.g
reads_file = args.r

check_input(genome_file, reads_file)
indexing(genome_file, "indexing")
mapping("indexing", reads_file, "mapping.sam")
sam_to_bam("mapping.sam", "mapping.bam") 
assembling("mapping.bam", "assembly.gtf")
orfsearching("assembly.gtf", genome_file)
