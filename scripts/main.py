#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import yaml

#VARIABLES 
cwd = os.getcwd()

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

#dont know if needed 
def fastq_to_fasta(fastqFile):
    output_file = fastqFile.split(".")[0] + ".fasta"
    try:
        seqtk_command = [
            "/home/s-amknut/GALBA/tools/seqtk/seqtk", #Neuen Pfad angeben
            "seq",
            "-a",
            fastqFile,
            "> ",
            output_file
        ]

        print("Converting .fastq to .fasta file...")

        with open(output_file, 'w') as output:
            result = subprocess.run(seqtk_command, stdout=output, capture_output=False)

        if result.returncode == 0:
            print("Conversion from .fastq to .fasta file completed successfully")
            return output_file

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

def file_name(path):
    path = path.split("/")
    name_with_dot = path[-1]
    name = name_with_dot.split(".")
    if name[0].endswith("_1"):
        name = name[0].split("_1")
        return name[0]
    else:
        return name[0]

def file_name_paired(path):  #not used yet
    path = path.split("/")
    name_with_dot = path[-1]
    name_with_number = name_with_dot.split(".")
    if name_with_number.endswith("_1"):
        name = name_with_number.split("_1")
        return name_with_number[0]
    elif name_with_number.endswith("_2"):
        name = name_with_number.split("_2")
        return name_with_number[0]
    else:
        return name[0]

def file_name_and_format(path): #not used yet
    path = path.split("/")
    return path[-1]

def first_or_second(path):
    path = path.split("/")
    name_with_dot = path[-1]
    name = name_with_dot.split(".")
    if name[0].endswith("1"):
        return "first"
    if name[0].endswith("2"):
        return "second"

def indexing(genome_fasta, output):
    hisat2_build_path = "/opt/hisat2/hisat2-build"
    #hisat2_build_path = "home/s-amknut/GALBA/tools/hisat2/hisat2-build"
    #output = file_name(genome_fasta) 
    try:
        hisat2_build_command = [
            hisat2_build_path, 
            "--quiet",
            "-p",
            "4",        
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

def mapping_short(genome, rnaseq_single_sets, rnaseq_paired_sets):
    #genome = file_name(genome)
    hisat2_path = "/opt/hisat2/hisat2"
    if file_format(rnaseq_single_sets[0]) == "fasta": #Muss ich hier beide sets testen?
        format_option = "-f"
    if file_format(rnaseq_single_sets[0]) == "fastq":
        format_option = ""
    try:
        for s in rnaseq_single_sets: 
            set_name = file_name(s) #-1 und -2 entfernen
            output_sam = set_name + ".sam" #Samfile noch in current working directory speichern
            hisat2_command = [
                    hisat2_path, 
                    format_option,                       
                    "-x", genome,            
                    "-U", s,
                    "--dta",
                    "-p", str(threads),
                    "-S", output_sam,                  
                ]
            print("Mapping set " + file_name(s) + " to genome...")
            
            result = subprocess.run(hisat2_command, capture_output=True)
            
            if result.returncode == 0:
                print("Mapping of set " + file_name(s) + " completed successfully")
            else:
                print("Error during mapping of set " + file_name(s) + ":")
                print(result.stderr)

        for s in rnaseq_paired_sets: 
            if first_or_second(s) == "first":
                set1 = s
                output_sam = file_name(s) + ".sam"
                continue
            if first_or_second(s) == "second":
                set2 = s

            hisat2_command = [
                    hisat2_path, 
                    format_option,                       
                    "-x", genome,            
                    "-1", set1,
                    "-2", set2,   
                    "--dta",
                    "-p", str(threads),
                    "-S", output_sam,                  
                ]
            print("Mapping set " + file_name(s) + " to genome...")
            
            result = subprocess.run(hisat2_command, capture_output=True)
            
            if result.returncode == 0:
                print("Mapping of set " + file_name(s) + " completed successfully")
            else:
                print("Error during mapping of set " + file_name(s) + ":")
                print(result.stderr)
            
    except Exception as e:
        print("Could not run hisat2 command for set: " + file_name(s)) 
        
    #if file_format(reads_file) == "fastq":
     #   fastq_to_fasta(reads_file)
      #  reads_file = "reads.fa"
    #try:
        #which Pfad finden
     #   hisat2_command = [
      #      hisat2_path,
       #     "-f",                        
        #    "-x", indexed_genome,            
         #   "-U", reads_file,            
          #  "-S", output_sam,            
           # "--no-spliced-alignment"        
        #]
        
        #print("Mapping reads to genome...")
        
        #result = subprocess.run(hisat2_command, capture_output=True)
        
        #if result.returncode == 0:
        #    print("Mapping completed successfully")
        #else:
         #   print("Error during mapping:")
          #  print(result.stderr)
    
    #except Exception as e:
     #   print("Could not run hisat2 command.") 

def mapping_long(genome, reads_long):
    try :
        for s in reads_long:
            set_name = file_name(s)
            output_sam = set_name + ".sam"
            minimap2_command = [
                "opt/minimap2/minimap2", 
                "-a", #NICHT DOCH RICHTIG?
                genome,
                s,
                "-o",
                output_sam
            ]
    
            print("Mapping isoseq set: " + s +" to genome...")
            #print("Running command:", "".join(minimap2_command))

            result = subprocess.run(minimap2_command, capture_output=True)

            if result.returncode == 0:
                print("Long read mapping of set " + set_name + " completed successfully")
            else:
                print("Error during long read mapping:")
                print(result.stderr)
    
    except Exception:
        print("Could not run minimap2 command.")

def sam_to_bam(rna_paired_sets, rna_single_sets, rna_long_sets):
    try:    
        #combined_lists = rna_paired_sets[0::2] + rna_single_sets + rna_long_sets 
        combined_lists = rna_paired_sets[0::2]
        bam_file_list = []
        for s in combined_lists:
            set_name = file_name(s)
            sam_file = set_name + ".sam"
            output_bam = set_name + ".bam"
            samtools_command = [
                "samtools",
                "sort",
                sam_file,
                "-o",
                output_bam
            ]

            print("Converting "+ set_name + ".sam to .bam file...")
            #result = subprocess.run(samtools_command, capture_output=True)

            #if result.returncode == 0:
            print("Conversion from .sam to .bam file completed successfully")
            bam_file_list.append(output_bam)

            #else:
             #   print("Error during conversion:")
              #  print(result.stderr)
        print (bam_file_list)
        return bam_file_list

    except Exception:
        print("Could not run samtools command.")

#Input überarbeiten, sodass short und long einzeln eingegeben werden kann
def mergeBamFiles(bam_file_list):
    bam_string = ""
    for file in bam_file_list:
        bam_string = bam_string + " " + file #Ein Leerzeichen am Anfang zu viel
    print(bam_string)
    try:
        print("Merging bam files...")
        command = "samtools merge -f -o entireMapping.bam entireRna896_1.bam entireRna664_1.bam"
        command1 = [
            "samtools",
            "merge",
            "-f"
            "-o",
            "entireMapping.bam",
            "entireRna896_1.bam",
            "entireRna664_1.bam"
        ]
        result = os.system(command) #Reminder: Didnt work with subprocess.run, maybe need a better solution?
        #result = subprocess.run(command1, capture_output=True) 

        if result== 0:
            print("Merged bam files successfully")

        else:
            print("Error during merging bam files")

    except Exception:
        print("Could not run samtools command.")

def assembling(shortBamFile, longBamFile, output_gtf):
    try:
        print("Assembling the reads...")
        #command = "stringtie -p "+str(threads)+" -o "+output_gtf+" "+bamFile
        command = "stringtie --mix -o " + output_gtf + " " + shortBamFile + " " + longBamFile
        #command = "stringtie " + shortBamFile + " " + longBamFile + " -o " + output_gtf
        print("Running command:", command)
        result = os.system(command) #Reminder: Didnt work with subprocess.run, maybe need a better solution?

      #  if result== 0:
       #     print("Assembled reads successfully")

       # else:
        #    print("Error during Assembly")

    except Exception:
        print("Could not run stringtie command.")

def assembling(rnaseq_paired_sets, rnaseq_single_sets, isoseq_sets):
    try:
        print("Assembling the reads...")
        combined_lists = rnaseq_paired_sets + rnaseq_single_sets
        for s in combined_lists:
            set_name = file_name(s)
            bam_file = set_name + ".bam"
            output_gtf = set_name + ".gtf"
            command = "stringtie -p "+str(threads)+" -o " + output_gtf + " " + bam_file

            result = os.system(command) #Reminder: Didnt work with subprocess.run, maybe need a better solution?

            if result == 0:
                gtf_list = gtf_list.append(output_gtf)
                print("Assembled reads successfully")

            else:   
                print("Error during Assembly")

    except Exception:
        print("Could not run stringtie command for short reads.")
    
    try:
        for s in isoseq_sets:
            set_name = file_name(s)
            bam_file = set_name + ".bam"
            output_gtf = set_name + ".gtf"
            command = "stringtie -L -p "+str(threads)+" -o " + output_gtf + " " + bam_file
        
            print("Running command:", command)
            result = os.system(command) #Reminder: Didnt work with subprocess.run, maybe need a better solution?

            if result== 0:
                print("Assembled reads successfully")
                gtf_list = gtf_list.append(output_gtf)

            else:
                print("Error during Assembly")
        
        print("Merge Transcripts...")
        command = "stringtie --merge -o transcripts_merged.gtf" +  gtf_list
        result = os.system(command)
        if result== 0:
            print("Merged transcripts successfully")
        else:
            print("Error during merging transcripts")

    except Exception:
        print("Could not run stringtie command for isoseq reads.")


def orfsearching(assembly_gtf, genome_fa, output_fa):
    try:
        gffread_command = [
            "/opt/gffread/gffread",
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
            "/opt/eviann/TransDecoder.Predict",
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

def load_config(config_file):
    with open(config_file, "r") as config_file:
        input_files = yaml.safe_load(config_file)
        return input_files

#MAIN
parser = argparse.ArgumentParser()  
parser.add_argument('-t', default=4, help='Number of threads', required=False)
parser.add_argument('-y', help='Config file input', required=True)

args = parser.parse_args()
#genome_file = args.g #50.000 von 1.985.779
#reads_short = args.s #100.000 von 87.429.668
threads = args.t

input_files = load_config(args.y)
genome_file = input_files["genome"]
rnaseq_paired_sets = input_files["rnaseq_paired_sets"]
rnaseq_single_sets = input_files["rnaseq_single_sets"]
isoseq_sets = input_files["isoseq_sets"]

#check_input(genome_file, reads_file) hier nochmal gut Lösung überlegen
indexing(genome_file, "genome")
print("mapping_short bis mergeBamFiles")
#mapping_short(genome_file, rnaseq_single_sets, rnaseq_paired_sets)
#mapping_long(genome_file, isoseq_sets)  
#bam_file_list = sam_to_bam(rnaseq_paired_sets, rnaseq_single_sets, isoseq_sets) 
#mergeBamFiles(bam_file_list)
#assembling(rnaseq_paired_sets, rnaseq_single_sets, isoseq_sets) #transcripts_merged.gtf not here!!
#orfsearching("transcripts_merged.gtf", genome_file, "transcripts.fasta")

#TO DOs:
#-submit.sh -B ändern
#-Variablen und Funktionsnamen anpassen
#-Nur input[isoseq] und co wenn diese "Kategorie" auch in der config file vorhanden ist
#-Funktion die prüft ob files vorhanden wie CreateThis() von GeneMark 
#-CheckInput Funktion anpassen und logisch machen 
#-cwd integrieren
#-verstehen was das --dta in hisat2 bedeutet 
#-prints überarbeiten
#--Ausgaben für subprocess.run überarbeiten, bei error mehr ausprinten lassen.

#FRAGEN:
#-Vor jedem Aufruf alte files löschen?
#-Genome Namen bei indexing beibehalten oder "genome" nennen?
#-Minimap2 option nicht doch richtig?
#-Ablauf richtig: 2x mergen aufrufen für einmal short und einmal long und dann stringtie mit --mix aufrufen?
