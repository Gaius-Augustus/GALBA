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
        print("Converting reads file to fasta format...")
        fastq_to_fasta(reads_file)
    if file_format(genome_file) == "fastq":
        print("Input genome file is in fastq format.")
        print("Converting genome file to fasta format...")
        fastq_to_fasta(genome_file)
    if file_format(reads_file) == "unknown":
        sys.exit("Error Reads file: Unknown file format")
    if file_format(genome_file) == "unknown":
        sys.exit("Error Genome file: Unknown file format")

def fastq_to_fasta(fastqFile):
    try:
        seqtk_command = [
            "/home/s-amknut/GALBA/tools/seqtk/seqtk",
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
            sys.exit(1)
    
    except Exception:
        print("Could not run seqtk command.")
        sys.exit(1)

def file_format(file):
    with open(file, 'r') as f:
        #read and save first line in string & remove whitespaces with strip()
        #every next call of .readline() will read the next line
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
    hisat2_build_path = "/opt/hisat2/hisat2-build"
    #hisat2_build_path = "/home/s-amknut/GALBA/tools/hisat2/hisat2-build"
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

def mapping_short(output_indexing, reads_file, output_sam):
    hisat2_path = "/opt/hisat2/hisat2"
    #hisat2_path = "/home/s-amknut/GALBA/tools/hisat2/hisat2"
    if file_format(reads_file) == "fastq":
        fastq_to_fasta(reads_file)
        reads_file = "reads.fa"
    try:
        #which Pfad finden
        hisat2_command = [
            hisat2_path,
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
        print("Could not run hisat2 command.") 

def mapping_long(genome, reads_long, output_sam):
    try :
        minimap2_command = [
            "/home/s-amknut/GALBA/tools/minimap2/minimap2",
            "-a",
            genome,
            reads_long,
            ">",
            output_sam
        ]
    
        print("Mapping long reads to genome...")
        print("Running command:", "".join(minimap2_command))

        result = subprocess.run(minimap2_command, capture_output=True)

        if result.returncode == 0:
            print("Long read mapping completed successfully")
        else:
            print("Error during long read mapping:")
            print(result.stderr)
    
    except Exception:
        print("Could not run minimap2 command.")

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
        command = "stringtie -p "+str(cores)+" -o "+output_gtf+" "+bamFile
        print("Running command:", command)
        result = os.system(command) #Reminder: Didnt work with subprocess.run, maybe need a better solution?

        if result== 0:
            print("Assembled reads successfully")

        else:
            print("Error during Assembly")

    except Exception:
        print("Could not run stringtie command.")


def orfsearching(assembly_gtf, genome_fa, output_fa):
    try:
        gffread_command = [
            "/home/s-amknut/GALBA/tools/gffread/gffread",
            "-w",
            output_fa,
            "-g",
            genome_fa,
            assembly_gtf
        ]
        print("Preparing Transcripts using gffread...")
        #print("Running command:", "".join(gffread_command))
        result = subprocess.run(gffread_command, capture_output=True)

        if result.returncode == 0:
            print("Transcripts prepared successfully")
        else:
            print("Error during preparing transcripts")
            print(result.stderr)    
    
    except Exception:
        print("Could not run gffread command.")
        sys.exit(1)

    try:
        trans = "transcripts.fasta"
        path = "/home/s-amknut/GALBA/tools/eviann/TransDecoder.LongOrfs"
        print("Searching for ORFs...")
        longORF_command = [
            "TransDecoder.LongOrfs",
            "-t",
            trans
        ]
       
        #print("Running command:", "".join(transdecoder_command))
        result= subprocess.run(longORF_command, capture_output=True)
        if result.returncode == 0:
            print("ORFs found successfully")
        else:
            print("Error during ORF search")
            print(result.stderr)
    except Exception:
        print("Could not run TransDecoder command.")

    try:
        trans = "transcripts.fasta"
        print("Searching for ORFs...")
        predict_command = [
            "/home/s-amknut/GALBA/tools/eviann/TransDecoder.Predict",
            "-t",
            trans
        ]
       
        #print("Running command:", "".join(transdecoder_command))
        result= subprocess.run(predict_command, capture_output=True)
        if result.returncode == 0:
            print("Predict successfully")
        else:
            print("Error during ORF search")
            print(result.stderr)
    except Exception:
        print("Could not run TransDecoder command.") 


#MAIN
parser = argparse.ArgumentParser()  
parser.add_argument('-g', help='Genome file', required=True)
parser.add_argument('-s', help='Short reads file', required=False) #entweder/oder programmieren
parser.add_argument('-l', help='Long reads file', required=False)
parser.add_argument('-t', help='Number of threads', required=False)

args = parser.parse_args()
genome_file = args.g #50.000 von 1.985.779
reads_short = args.s #100.000 von 87.429.668
reads_long = args.l
threads = args.t

#check_input(genome_file, reads_file) hier nochmal gut Lösung überlegen
indexing(genome_file, "genome")
mapping_short("genome", reads_short, "mapping_short.sam")
mapping_long(genome_file, reads_long, "mapping_long.sam")
#sam_to_bam("mapping.sam", "mapping.bam") 
#assembling("mapping.bam", "assembly.gtf")
#orfsearching("assembly.gtf", genome_file, "transcripts.fasta")


'''
def install_cpan():
    if not os.path.isfile("/usr/bin/cpan"):
        try:
            result = os.system("echo yes | cpan")
            if result.returncode == 0:
               print("CPAN configuration set automatically")
            else:
               print("Error during CPAN configuration")
               print(result.stderr)
        except Exception:
            print("Error with setting cpan configuration automatically.")
            sys.exit(1)

        try:
            print("Installing URI::Escape...")
            #result = subprocess.run(["install", "URI::Escape"], capture_output=True)
            result = os.system("cpan URI::Escape")
            result = os.system("sudo cpanm install URI::Escape")
            if result.returncode == 0:
                print("URI::Escape installed successfully")
                #subprocess.run(["exit"], capture_output=True)
            else:
                print("Error during installation of URI::Escape")
                print(result.stderr)

        #result = os.system("perl -MCPAN -e shell")

        #if result.returncode == 0:
         #   print("CPAN interactive shell opened")
        #else:
         #   print("Error during opening of CPAN interactive shell")
          #  print(result.stderr)
       # print("Check if capnminus is installed...")
        #result = os.system("which cpanm")
        #if result.returncode != 0:
         #   print("cpanm is not installed")
          #  print("Installing cpanminus...")
           # result = os.system("curl -L https://cpanmin.us | perl - --sudo App::cpanminus")
            #if result.returncode != 0:
             #   print("Error during installation of cpanminus")
              #  print(result.stderr)
            #else:
             #   print("cpanminus installed successfully")
        #result = subprocess.run(["perl", "-MURI::Escape", "-e", "print 'URI::Escape is installed'"], capture_output=True)
        #result = subprocess.run(["perl", "-MCPAN", "-e", "shell"], capture_output=True)
        #result = subprocess.run(["cpan", "App::cpanminus"], capture_output=True)
        #result = os.system("cpan App::cpanminus")
        except Exception:
            print("Could not install cpan.")
            sys.exit(1)
    
def install_URI_escape():
    try:
        if not os.path.isfile("/usr/bin/cpanm"):
            print("cpanm is not installed")
            print("Installing cpanminus...")
            result = os.system("apt install cpanminus")
            if result.returncode != 0:
                print("Error during installation of cpanminus")
                print(result.stderr)
            else:
                print("cpanminus installed successfully")
        else:
            print("cpanminus is already installed")
            #if not
        #result = os.system("echo yes | cpan")
        #if result.returncode == 0:
         #   print("CPAN configuration set automatically")
        #else:
         #   print("Error during CPAN configuration")
          #  print(result.stderr)

        #result = os.system("perl -MCPAN -e shell")

        #if result.returncode == 0:
         #   print("CPAN interactive shell opened")
        #else:
         #   print("Error during opening of CPAN interactive shell")
          #  print(result.stderr)
        
    
    
        print("Installing URI::Escape...")
        #result = subprocess.run(["install", "URI::Escape"], capture_output=True)
        result = os.system("cpanm URI::Escape")
        if result.returncode == 0:
            print("URI::Escape installed successfully")
            #subprocess.run(["exit"], capture_output=True)
        else:
            print("Error during installation of URI::Escape")
            print(result.stderr)
    
    except Exception:
        print("Could not install URI::Escape.")
        sys.exit(1)
'''