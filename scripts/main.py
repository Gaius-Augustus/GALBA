#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import re
import yaml
import pandas as pd
import math
from Bio import SeqIO

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

def file_suffix(path):  #not used yet
    file = path.split("/")[-1]
    return file.split(".")[-1]

def file_name(path):
    file = path.split("/")[-1]
    name = file.split(".")[0]
    if name.endswith("_1") or name.endswith("_2"):
        prefix = name.split("_")[0]
        return prefix
    else:
        return name

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

def indexing(genome_fasta):
    hisat2_build_path = "hisat2-build"
    try:
        hisat2_build_command = [
            hisat2_build_path, 
            "--quiet",
            "-p",
            str(threads),        
            genome_fasta,           
            "genome"     
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

def mapping_short(rnaseq_paired_sets, rnaseq_single_sets):
    alignments_list = []
    try:
        if not rnaseq_single_sets == []:
            if file_format(rnaseq_single_sets[0]) == "fasta": #Muss ich hier beide sets testen?
                format_option = "-f"
            if file_format(rnaseq_single_sets[0]) == "fastq":
                format_option = ""
            string_with_sets = ",".join(rnaseq_single_sets)
            output_1 = "alignment_single_rnaseq_test1.sam" 
            hisat2_command = [
                    "hisat2",  
                    format_option,                       
                    "-x", "genome",            
                    "-U", string_with_sets,
                    "--dta",
                    "-p", str(threads),
                    "-S", output_1,                  
                ]
            print("Mapping single-end RNAseq to genome...")
            
            result = subprocess.run(hisat2_command, capture_output=True)
            
            if result.returncode == 0:
                print("Mapping of set of single-end RNASeq evidence completed successfully")
                alignments_list.append(output_1)

            else:
                print("Error during mapping of single-end rnaseq data:")
                print(result.stderr)
        
    except Exception as e:
        print("Could not run hisat2 command for single-end RNA-seq data.") 
        exit(1)

    try:
        if not rnaseq_paired_sets == []:
            if file_format(rnaseq_paired_sets[0]) == "fasta":
                format_option = "-f"
            if file_format(rnaseq_paired_sets[0]) == "fastq":
                format_option = ""
            string_with_first = ",".join(rnaseq_paired_sets[0::2])
            string_with_second = ",".join(rnaseq_paired_sets[1::2])
            output_2 = "alignment_paired_rnaseq_test1.sam" 
            hisat2_command = [
                    "hisat2", 
                    format_option,                       
                    "-x", "genome",            
                    "-1", string_with_first,
                    "-2", string_with_second,   
                    "--dta",
                    "-p", str(threads),
                    "-S", output_2,                  
                ]
            print("Mapping paired-end rnaseq data to genome...")
            
            result = subprocess.run(hisat2_command, capture_output=True)
            
            if result.returncode == 0:
                print("Mapping of paired-end rnaseq data completed successfully")
                alignments_list.append(output_2)
            else:
                print("Error during mapping of paired-end rnaseq data:")
                print(result.stderr)

    except Exception:
        print("Could not run hisat2 command for paired-end RNA-seq data.")
        exit(1)
    return alignments_list 

def mapping_long(genome, isoseq_sets):
    try :
        output_sam = "alignment_isoseq_test1.sam" #threads neu prüfen
        minimap2_command = ["minimap2", "-ax", "splice", "-uf", "-C5", genome, "-t", str(threads)] + isoseq_sets + ["-o", output_sam]
        #Threads noch hinzufügen
        #We can use -C5 for reads with low error rates like isoseq 

        print("Mapping isoseq sets to genome...")
        print("Running command:", " ".join(minimap2_command))

        result = subprocess.run(minimap2_command, capture_output=True)

        if result.returncode == 0:
            print("Mapping of isoseq reads completed successfully")
            return output_sam
        else:
            print("Error during mapping of isoseq data:")
            print(result.stderr)

    except Exception:
        print("Could not run minimap2 command.")

def sam_to_bam(sam_file_list):
    try:    
        for samfile in sam_file_list:
            output_bam = file_name(samfile) + ".bam"
            samtools_command = [
                "samtools",
                "sort",
                "-@",
                str(threads), #NEU THREADS HINZUGEFÜGT
                samfile,
                "-o",
                output_bam
            ]

            print("Converting " + samfile +" to " + output_bam + "...")
            result = subprocess.run(samtools_command, capture_output=True)

            if result.returncode == 0:
                print("Conversion from .sam to .bam file completed successfully")

            else:
                print("Error during conversion:")
                print(result.stderr)

    except Exception:
        print("Could not run samtools command.")

def merge_bam_files(bamfile_1, bamfile_2): 
    try:
        print("Merging bam files " + bamfile_1 + " and " + bamfile_2 + "...")
        output_bam = "alignment_merged_rnaseq_test1.bam"
        command1 = [
            "samtools",
            "merge",
            "-f", #-f is to override the output file if it already exists
            "-o",
            output_bam,
            bamfile_1,
            bamfile_2
        ] #keine threads, denn sequenziell nicht parallel
        result = subprocess.run(command1, capture_output=True) 
        #print(result.stdout)
        #print(result.stderr)
        if result.returncode== 0:
            print("Merged rnaseq bam files successfully")
            return output_bam
        else:
            print("Error during merging of rnaseq bam files")

    except Exception:
        print("Could not run samtools command.")

def assembling(alignment_rnaseq, alignment_isoseq):
    try:
        print("Assembling the reads...")
        output_gtf = "transcripts_mixed_test1.gtf"
        #alignment_rnaseq = "alignment_paired_rnaseq_test1.bam"
        #alignment_isoseq = "alignment_isoseq.bam"

        if args.mixed:
            command_mixed = [  #-p str(threads) noch hinzufügen
                "stringtie",
                "-p",
                str(threads),
                "-o",
                output_gtf,
                "--mix",
                alignment_rnaseq,
                alignment_isoseq
            ]

            result = subprocess.run(command_mixed, capture_output=True)

            if result.returncode == 0:
                print("Assembled rnaseq and isoseq reads successfully")

            else:
                print("Error during Assembly of rnaseq and isoseq reads")

        if args.rnaseq:
            print("In args.rnaseq von assmebling")
            command_rnaseq = [
                "stringtie",
                "-p",
                str(threads),
                "-o",
                output_gtf,
                alignment_rnaseq
            ]
            result = subprocess.run(command_rnaseq, capture_output=True)

            if result.returncode == 0:
                print("Assembled rnaseq reads successfully")
            else:
                print("Error during Assembly of rnaseq reads")
            
        if args.isoseq:
            command_isoseq = [
                "stringtie",
                "-L",
                "-p",
                str(threads),
                "-o",
                output_gtf,
                alignment_isoseq
            ]

            result = subprocess.run(command_isoseq, capture_output=True)

            if result.returncode == 0:
                print("Assembled isoseq reads successfully")
            else:
                print("Error during Assembly of isoseq reads")

    except Exception:
        print("Could not run stringtie command.")


def orfsearching(genome_fa, transcripts_gtf):
    try:
        output_fa = "transcripts_test1.fasta" 
        gffread_command = [
            "gffread",
            "-w",
            output_fa,
            "-g",
            genome_fa,
            transcripts_gtf
        ]
        print("Construct a fasta file from the transcripts.gtf file given by Stringtie...")
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

    #No threads option for Transdecoder (ChatGPT says option would be to split input files)
    try:
        print("Extract the long open reading frames...")
        longORF_command = [
            "TransDecoder.LongOrfs",
            "-t",
            output_fa
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

    #Optional on Transdecoder paige: Predicting ORFs with homology to known proteins

    try:
        print("Predict the likely coding regions...")
        predict_command = [
            "TransDecoder.Predict",
            "-t",
            output_fa
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

def protein_aligning(genome, protein, alignment_scoring):
    try: 
        #output_aln = "protein_alignment.aln"
        command = [
            "miniprot",
            "-t",
            str(threads),
            genome,
            protein,
            "--aln"
            "-o",
            "miniprot.aln"
        ]
        print("Aligning proteins to genome...")

        result = subprocess.run(command, capture_output=True)
        if result.returncode == 0:
            print("Proteins aligned successfully")
        else:
            print("Error during protein alignment with miniprot")
            print("stdout:", result.stdout.decode())
            print("stderr:", result.stderr.decode())

    except Exception:
        print("Could not run miniprot command.")
        sys.exit(1)

    try: 
        command = f"miniprot_boundary_scorer -o miniprot_parsed.gff -s {alignment_scoring} < miniprot.aln"
        print("Running command:", command)
        print("Scoring the alignment...")
        
        result = subprocess.run(command, shell=True, capture_output=True)

        if result.returncode == 0:
            print("Alignment scored successfully with miniprot_boundary_scorer!")
        else:
            print("Error during scoring the alignment with miniport_boundary_scorer!")
            print("stdout:", result.stdout.decode())
            print("stderr:", result.stderr.decode())
    
    except Exception:
        print("Could not run miniprot_boundary_scorer command.")
        sys.exit(1)

    try:
        command = [
            "miniprothint.py",
            "miniprot_parsed.gff",
            "--workdir",
            "miniprothint"
        ]

        print("Creating hints for Augustus...")

        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Hints created successfully")
        else:
            print("Error during creating hints by miniprothint")
            print("stdout:", result.stdout.decode())
            print("stderr:", result.stderr.decode())
            print(result.stderr)
        
    except Exception:
        print("Could not run miniprothint command.")
        sys.exit(1)

def shorten_incomplete_Orfs(transdecoder_pep):
    try:
        short_start_dict = {}
        with open("shortened_candidates.pep", "w") as output:
            for record in SeqIO.parse(transdecoder_pep, "fasta"):
                if "type:5prime_partial" in record.description or "type:internal" in record.description:
                    m_position = record.seq.find("M")
                    if m_position != -1:
                        record.seq = record.seq[m_position:]
                        description = record.description.split(" ")
                        coords = re.search(r":(\d+)-(\d+)\([\+\-]\)", description[7])
                        new_length = len(record.seq)
                        if coords:
                            old_start = int(coords.group(1)) 
                            new_start = old_start + m_position #Nicht -1, denn m_position ist 0-based
                            stop = int(coords.group(2))
                            description = re.sub(f"{old_start}-{stop}", f"{new_start}-{stop}", record.description)
                            description = re.sub(r"len:\d+", f"len:{new_length}", description)
                            record.description = description
                            SeqIO.write(record, output, "fasta")
                            short_start_dict[record.id] = m_position
        return short_start_dict
    
    except Exception:
        print("Error during writing shortened candidates .pep file.")
        sys.exit(1)

def make_diamond_db(protein_file):
    try:
        command = [
            "diamond",
            "makedb",
            "--in",
            protein_file,
            "-d",
            "protein_db"
        ]
        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Database created successfully")
        else:
            print("Error during creating database")
            print(result.stderr)
    
    except Exception:
        print("Could not run diamond makedb command.")
        sys.exit(1)

def validating_ORFs(transdecoder_pep, output_tsv):
    try:
        command = [
            "diamond",
            "blastp",
            "-d",
            "protein_db",
            "-q",
            transdecoder_pep,
            "-o",
            output_tsv
        ]
        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Blastp search for completed successfully and ", output_tsv, " was created.")
        else:
            print("Error during diamond blastp search.")
            print(result.stderr)
    
    except Exception:
        print("Could not run diamond blastp command for complete CDS candidates.")
        sys.exit(1)

def get_cds_classification(normal_tsv, shortened_tsv, short_start_dict):
    #Weil hier nach CDS-ID und nach Protein gemerged wird, werden einzelne herausgefiltert, die nicht in beiden vorkommen. 
    header_list = ["cdsID", "proteinID", "percIdentMatches", "alignLength", "mismatches", "gapOpenings", "alignStart", "alignEnd", "proteinStart", "proteinEnd", "eValue", "bitScore"]
    df_shortened = pd.read_csv(shortened_tsv, delimiter='\t', header=None, names=header_list)
    df_normal = pd.read_csv(normal_tsv, delimiter='\t', header=None, names=header_list)
    merged_df = pd.merge(df_shortened, df_normal, on=["cdsID", "proteinID"], suffixes=('_short', '_normal'))
    merged_df = merged_df.drop(columns=["alignLength_short", "alignLength_normal", "mismatches_short", "mismatches_normal", "gapOpenings_short", "gapOpenings_normal", "alignEnd_short", "alignEnd_normal", "proteinEnd_short", "proteinEnd_normal", "eValue_short", "eValue_normal"])
    merged_df["shortStart"] = merged_df["cdsID"].map(short_start_dict)
    merged_df["supportScore"] = None
    #merged_df["classification"] = None

    for i, cds in merged_df.iterrows():
        q_incomplete_start = cds["proteinStart_normal"]
        t_incomplete_start = cds["alignStart_normal"]
        t_complete_start = cds["alignStart_short"] + cds["shortStart"]  #not -1 because alignmentstart is 1-based but Mposition not. +shortstart, weil Differenz von normal zu short gebraucht wird
        aai_incomplete = cds["percIdentMatches_normal"]
        aai_complete = cds["percIdentMatches_short"]
        if aai_complete == 0:
            aai_complete = 0.0001
        if aai_incomplete == 0:
            aai_incomplete = 0.0001
        match_log = math.log(aai_incomplete/aai_complete)
        #maybe betrag nehmen
        support_score = (t_complete_start - t_incomplete_start) - (q_incomplete_start - 1) + match_log**1000
        merged_df.at[i, "supportScore"] = support_score
        #if t_complete_start < t_incomplete_start:
            #print("cdsID: ", cds["cdsID"] , "incomplete Start: ", t_incomplete_start, "complete Start: ", t_complete_start)
    
    merged_df["bitScore_max"] = merged_df[["bitScore_short", "bitScore_normal"]].max(axis=1) #FRAGE: Hier richtig über bitscore das beste alignment zu finden?
    grouped = merged_df.groupby("cdsID")
    #from here on it takes a lot of runtime because of .append() and .iterrows()
    classifications = {}
    for groupname, groupdata in grouped:
        sorted_group = groupdata.sort_values(by="bitScore_max", ascending=False)
        incomplete = False
        for i, cds in sorted_group.head(25).iterrows():
            if cds["supportScore"] > 0:
                incomplete = True
                break
        if incomplete:
            classifications[cds["cdsID"]] = "incomplete" #112 13.8%
        else:
            classifications[cds["cdsID"]] = "complete" #699 86.2%
            #print("t_complete_start: ", cds["alignStart_normal"], "t_incomplete_start: ", cds["alignStart_short"] )

    pd.set_option('display.max_columns', None)
    #print(merged_df.head())
    #print(classifications)
    #print(count)
    return classifications

def get_hc_cds(diamond_tsv, classifications):
    header_list = ["cdsID", "proteinID", "percIdentMatches", "alignLength", "mismatches", "gapOpenings", "alignStart", "alignEnd", "proteinStart", "proteinEnd", "eValue", "bitScore"]
    df = pd.read_csv(diamond_tsv, delimiter='\t', header=None, names=header_list)
    q_length_dict = {}
    for record in SeqIO.parse(protein_file, "fasta"):
        q_length_dict[record.id] = len(record.seq)
    t_length_dict = {}
    for record in SeqIO.parse("revised_candidates.pep", "fasta"):
        t_length_dict[record.id] = len(record.seq)
    with open(diamond_tsv, "r") as tsv:
       # with open("hc_genes.pep", "w") as output:
        for line in tsv:    #7580 Elemente in tsv
            part = line.strip().split('\t')
            cdsID = part[0]
            #proteinID = part[1]
            #percIdentMatches = float(part[2])
            align_length = int(part[3])
            #mismatches = int(part[4])
            #gapOpenings = int(part[5])
            t_start = int(part[6])
            t_end = int(part[7])
            q_start = int(part[8])
            q_end = int(part[9])
            #eValue = float(part[10])
            #bitScore = float(part[11])
            q_length = q_length_dict[cdsID]
            t_length = t_length_dict[cdsID]
            
            start_condition = q_start - t_start < 6
            stop_condition = (q_length - q_end) - (t_length - t_end) < 21
            if classifications[cdsID] == "complete":
                if start_condition and stop_condition:
                    hc = True 

            if classifications[cdsID] == "3prime_partial" or classifications[cdsID] == "internal":
                continue
    
def get_optimized_pep_file(normal_pep, shortened_pep, classifications, short_start_dict):
    with open("revised_candidates.pep", "w") as output:
        shortened_pep_dict = {record.id: (record.seq, record.description) for record in SeqIO.parse(shortened_pep, "fasta")}
        classifications_for_hc = {}
        for record in SeqIO.parse(normal_pep, "fasta"):
            id = record.id
            if "type:complete" in record.description:
                classification = "complete"
            if "type:5prime_partial" in record.description:
                classification = "5prime_partial"
            if "type:internal" in record.description:
                classification = "internal"
            if "type:3prime_partial" in record.description:
                classification = "3prime_partial"

            if classification == "5prime_partial" or classification == "internal":
                if id in classifications:
                    if classifications[id] == "incomplete":
                        SeqIO.write(record, output, "fasta")  
                    else:
                        seq = shortened_pep_dict[id][0]
                        record.seq = seq
                        description = re.sub(r"len:\d+", f"len:{len(seq)-1}", record.description) #-1 damit * nicht gezählt wird
                        coords = re.search(r":(\d+)-(\d+)\([\+\-]\)", description)
                        start_normal = int(coords.group(1))
                        start_short = start_normal + short_start_dict[id] * 3 #Nicht -1, denn m_position ist 0-based
                        stop = int(coords.group(2))
                        strand = coords.group(0)[-2]
                        new_coords = f'{start_short}-{stop}({strand})'
                        description = re.sub(r'(\d+-\d+\([\+\-]\))', new_coords, description)
                        if "type:5prime_partial" in record.description:
                            description = description.replace("type:5prime_partial", "type:complete")
                            classification = "complete"
                        else:
                            description = description.replace("type:internal", "type:3prime_partial")
                            classification = "3prime_partial"
                        record.description = description   
                        SeqIO.write(record, output, "fasta")
                else:
                    SeqIO.write(record, output, "fasta") #If there is no classification, there was no protein evidence found for the CDS 
                    classifications[id] = "incomplete"   #(drüber nachdenken ob es sein kann, dass es nur für den normalen aber dafür für den kurzen evidenz gibt. Auch dann kommt cds in merged_df nicht vor)
            else:
                SeqIO.write(record, output, "fasta")

            classifications_for_hc[id] = classification

        return classifications_for_hc

    #4691 sind in classifications nicht drin, aber in normal_pep als 5prime/ internal
    #5502 sind in 5prime/internal in normal_pep
    #310 sind in normal_pep als 5prime_partial oder internal, aber nicht in short_pep_dict, weil sie kein M enthalten!!!
    #4682 sind in normal_pep als 5prime_partial oder internal, aber nicht in short_tsv
    #820 elemente sind in short_tsv
    #9 Elemente sind in short_tsv aber nicht in classifications -> Vermutung: Protein aligniert nur mit short oder nur mit normal und alignment taucht damit nicht in mergeddf auf

def parse_transdecoder_file(transdecoder_pep):
    transdecoder_id_dict = {}
    for record in SeqIO.parse(transdecoder_pep, "fasta"):
        id = record.id
        id = id.split(".p")[0]
        description = record.description.split(" ")
        cds_transcript_coords = re.search(r":(\d+)-(\d+)\((\+|\-)\)", description[7])
        #if cds_transcript_coords:
        start_cds_transcript = int(cds_transcript_coords.group(1))
        stop_cds_transcript = int(cds_transcript_coords.group(2))
        cds_length = int(stop_cds_transcript) - int(start_cds_transcript) + 1
        triple = (start_cds_transcript, stop_cds_transcript, cds_length)
        if id not in transdecoder_id_dict:
            transdecoder_id_dict[id] = [triple]
        else:
            transdecoder_id_dict[id].append(triple)
            
    filtered_dict = {id: val for id, val in transdecoder_id_dict.items() if len(val) == 1}
    return filtered_dict

def from_transcript_to_genome_coords(stringtie_gtf, transdecoder_id_dict): #better name: creating_annotation_file
    annotation_file = "annotation.gtf"
    with open(stringtie_gtf, "r") as stringtie, open(annotation_file, "w") as output:
        exon_coords_list = []
        for line in stringtie:
            if line.startswith("#"):
                continue
            else:                
                part = line.strip().split('\t')
                seqname = part[0]
                source = part[1]
                feature = part[2]
                start_genome = part[3]
                stop_genome = part[4]
                score = part[5]
                strand = part[6]
                frame = part[7]
                attributes = part[8]

                gene_id = re.search(r'gene_id "([^"]+)"', attributes)
                gene_id = gene_id.group(1)
                transcript_id = re.search(r'transcript_id "([^"]+)"', attributes)
                transcript_id = transcript_id.group(1)
                #Eventuell noch exon_number hinzufügen
                #Soll score mit rein? Habs erstmal rausgenommen
                
                if feature == 'transcript':
                    print("--------------Neues Transcript mit ID: ----------- ", transcript_id) 
                    output.write(f"{seqname}\tPreGalba\tgenome\t{start_genome}\t{stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    output.write(f"{seqname}\tPreGalba\ttranscript\t{start_genome}\t{stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    print("ANFANG: Transkriptlänge, cds Länge, cds index auf 0")
                    if len(exon_coords_list) > 1: #Für CDS muss es auch zugelassen sein, dass nur ein Exon vorhanden ist
                        transcript_id = exon_coords_list[0][2]
                        gene_id = exon_coords_list[0][3]
                        strand = exon_coords_list[0][4]
                        frame = exon_coords_list[0][5]
                        print("Exons für vorheriges Transkript vorhanden: ", transcript_id)
                        for i in range(len(exon_coords_list)-1):
                            curr_exon_stop = int(exon_coords_list[i][1])
                            next_exon_start = int(exon_coords_list[i+1][0])
                            intron_start = int(curr_exon_stop) + 1
                            intron_stop = int(next_exon_start) - 1
                            print("Intron von: ", intron_start, "bis: ", intron_stop)
                            output.write(f"{seqname}\tPreGalba\tintron\t{intron_start}\t{intron_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    if len(exon_coords_list) > 0:
                        if transcript_id in transdecoder_id_dict:
                            cds_transcript_start = int(transdecoder_id_dict[transcript_id][0][0])
                            cds_transcript_stop = int(transdecoder_id_dict[transcript_id][0][1])
                            cds_total_length = int(transdecoder_id_dict[transcript_id][0][2])
                            print("CDS ist insgesamt so lang: ", cds_total_length)
                            print("CDS Koordinaten in Transkript: ", cds_transcript_start, cds_transcript_stop)
                            curr_transcript_length = 0
                            cds_current_length = 0
                            for i in range(len(exon_coords_list)):
                                curr_exon_start = int(exon_coords_list[i][0]) 
                                curr_exon_stop = int(exon_coords_list[i][1])
                                print("Exonkoordinaten vom aktuellen Exon: ", curr_exon_start, curr_exon_stop)
                                if i > 0:
                                    curr_transcript_length += int(exon_coords_list[i-1][1]) - int(exon_coords_list[i-1][0]) + 1
                                    print("Vorherige Transkriptlänge wird hochgesetzt auf: ", curr_transcript_length)
                                #Anfang vom CDS:
                                if cds_current_length == 0:
                                    if curr_exon_start - curr_transcript_length + cds_transcript_start - 1 > curr_exon_stop: #GEÄNDERT
                                        x = curr_exon_start + cds_transcript_start - 1 #nochmal genau bestimmen GEÄNDERT
                                        print("CDS Startpunkt:", x , " liegt hinter Exonstoppunkt: ", curr_exon_stop, "Also zu nächstem Exon springen")
                                        fivePrimeUTR_start = curr_exon_start
                                        fivePrimeUTR_stop = curr_exon_stop
                                        if strand == "+":
                                            output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        else:
                                            output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        continue
                                    else: 
                                        cds_start_genome = curr_exon_start + int(cds_transcript_start) - curr_transcript_length - 1
                                        print("CDS Startpunkt liegt innerhalb des Exons bei: ", cds_start_genome)
                                        if cds_start_genome + cds_total_length - cds_current_length - 1 > curr_exon_stop: #GEÄNDERT
                                            print("CDS geht bis Exongrenze...")
                                            output.write(f"{seqname}\tPreGalba\tCDS\t{cds_start_genome}\t{curr_exon_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            print("CDS von: ", cds_start_genome, "bis: ", curr_exon_stop, "hinzugefügt")
                                            cds_current_length = curr_exon_stop - cds_start_genome + 1
                                            print("CDS aktuelle Länge: ", cds_current_length)
                                            fivePrimeUTR_start = curr_exon_start
                                            fivePrimeUTR_stop = cds_start_genome - 1
                                            start_codon_plus_start = cds_start_genome
                                            start_codon_plus_stop = cds_start_genome + 2
                                            if strand == "+":
                                                output.write(f"{seqname}\tPreGalba\tstart_codon\t{start_codon_plus_start}\t{start_codon_plus_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                if curr_exon_start != cds_start_genome:
                                                    output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            else:
                                                output.write(f"{seqname}\tPreGalba\tstop_codon\t{start_codon_plus_start}\t{start_codon_plus_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                if curr_exon_start != cds_start_genome:
                                                    output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")

                                        else:
                                            print("CDS innerhalb des Exons beendet...")
                                            cds_stop_genome = cds_start_genome + cds_total_length - 1
                                            cds_current_length = cds_total_length
                                            print("Aktueller Exonstoppunkt: ", stop_genome)
                                            print("CDS Stoppunkt: ", cds_stop_genome)
                                            output.write(f"{seqname}\tPreGalba\tCDS\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            fivePrimeUTR_start = curr_exon_start
                                            fivePrimeUTR_stop = cds_start_genome - 1
                                            threePrimeUTR_start = cds_stop_genome + 1
                                            threePrimeUTR_stop = curr_exon_stop
                                            print("CDS von: ", cds_start_genome, "bis: ", cds_stop_genome, "hinzugefügt")
                                            start_codon_plus_start = cds_start_genome
                                            start_codon_plus_stop = cds_start_genome + 2
                                            stop_codon_plus_start = cds_stop_genome - 2
                                            stop_codon_plus_stop = cds_stop_genome
                                            if strand == "+":
                                                if curr_exon_start != cds_start_genome:
                                                    output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                if curr_exon_stop != cds_stop_genome:
                                                    output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{threePrimeUTR_start}\t{threePrimeUTR_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                output.write(f"{seqname}\tPreGalba\tstart_codon\t{start_codon_plus_start}\t{start_codon_plus_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                output.write(f"{seqname}\tPreGalba\tstop_codon\t{stop_codon_plus_start}\t{stop_codon_plus_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            else:
                                                if curr_exon_start != cds_start_genome:
                                                    output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                if curr_exon_stop != cds_stop_genome:
                                                    output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{threePrimeUTR_start}\t{threePrimeUTR_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                output.write(f"{seqname}\tPreGalba\tstop_codon\t{start_codon_plus_start}\t{start_codon_plus_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                output.write(f"{seqname}\tPreGalba\tstart_codon\t{stop_codon_plus_start}\t{stop_codon_plus_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            continue
                                        
                                elif cds_current_length == cds_total_length:
                                    print("CDS bereits beendet. Rest der exons wird zu UTR")
                                    if strand == "+":
                                        output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{curr_exon_start}\t{curr_exon_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                    else:
                                        output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{curr_exon_start}\t{curr_exon_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                else:
                                    cds_start_genome = curr_exon_start
                                    if cds_current_length+int(curr_exon_stop)-int(curr_exon_start)+1<cds_total_length:
                                        print("CDS noch nicht beendet...")
                                        cds_stop_genome = int(curr_exon_stop) 
                                        print("CDS wird von Start des aktuellen Exons, bis Stopp des aktuellen Exons eingetragen")
                                        output.write(f"{seqname}\tPreGalba\tCDS\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        cds_current_length+=curr_exon_stop - curr_exon_start + 1
                                        print("Neue aktuelle Länge: ", cds_current_length)
                                    else:
                                        print("CDS innerhalb des Exons beendet...")
                                        cds_stop_genome = int(curr_exon_start) + cds_total_length - cds_current_length - 1 #GEÄNDERT (-1 hinzugefügt)
                                        print("Aktueller Exonstoppunkt: ", curr_exon_stop)
                                        print("CDS Stoppunkt: ", cds_stop_genome)
                                        output.write(f"{seqname}\tPreGalba\tCDS\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        threePrimeUTR_start = cds_stop_genome + 1
                                        threePrimeUTR_stop = curr_exon_stop
                                        stop_codon_plus_start = cds_stop_genome - 2
                                        stop_codon_plus_stop = cds_stop_genome
                                        if strand == "+":
                                            output.write(f"{seqname}\tPreGalba\tstop_codon\t{stop_codon_plus_start}\t{stop_codon_plus_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            if curr_exon_stop != cds_stop_genome:
                                                output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{threePrimeUTR_start}\t{threePrimeUTR_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        else:
                                            output.write(f"{seqname}\tPreGalba\tstart_codon\t{stop_codon_plus_start}\t{stop_codon_plus_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            if curr_exon_stop != cds_stop_genome:
                                                output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{threePrimeUTR_start}\t{threePrimeUTR_stop}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        continue
                                        
                    exon_coords_list.clear()
                        
                if feature == 'exon':
                    output.write(f"{seqname}\tPreGalba\texon\t{start_genome}\t{stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    exon_number = re.search(r'exon_number "([^"]+)"', attributes)
                    exon_number = exon_number.group(1)
                    exon_coords_list.append((start_genome, stop_genome, transcript_id, gene_id, strand, frame))
    
    '''
    command = [
        "sort",
        "-k9,9",
        "-t$'\t'",
        annotation_file,
        "-o",
        "annotation_sorted.gtf"
    ]
    subprocess.run(command)
'''

#(Zwischen cds im selben transkript liegen Utrs? Im selben Exon hört ein cds auf und anderes beginnt?)
#--> Schreibt man die mit auf?
#-Wie verwende ich die overlap.pl Datei? Ich muss ja eigentlich sagen können, wann zwei sets als gleich bezeichnet werden. key error
#-gffcompare versucht zu verwenden aber make funktioniert nicht
#-Wollen wir zusätzlich miniprot predictete ORFs verwenden? Stand jetzt nutzen wir miinniprot nur um zu schauen, dass Konflikte vermieden werden
#-Sind am Ende auch nur die HC Gene die die wir reinschreiben in die annot file oder schreiben wir alles rein?
                           
'''
Theorie:
-High confidence CDS finden mit folgenden Kriterien:
    Für complete:
        > Wenn CDS mit irgendeinem Protein die Gleichung erfüllt wird CDS als HC kategorisiert
    Für incomplete:
        > Wenn irgendein Protein C-terminus (alignment Ende) abdeckt (Wie in Gleichung), wird CDS als HC kategorisiert
        > Wenn das beste alignment nicht die Sequenz von Anfang an abdeckt, wird es gekürzt (?)
    Ohne Protein support oder ohne Gleichung zu erfüllen:
        > Betrachten nur die ohne Protein support
        > Nur das längste CDS wird pro Gen ausgewählt für dieses muss gelten: (Nur längste ohne Protein support?)
            > Länge >= 300 nucleotides 
            > GMS-T log-odds score > 50, #EIGENEN SCORE AUSWÄHLEN?
            > Complete
            > InFrame Stop Codon in the 5' UTR
            > No conflict with other gene prediction in same locus
            > No overlap with other gene prediction in same locus
                --> Mapping predicted gene to genome (?) and compare exon-intron structure with potentially conflicting MiniProthint prediction 
            
Plan:
-In Kürzungsfunktion instant_complete_dictonary erstellen (Für complete und 3' partial) 
-In Candidatesfunktion later_complete_dictonary mit gekürzte complete ORFS erstellen und incomplete_dictonary erstellen 
-In normal.tsv Gleichung für instant_complete_ORFS berechnen und HC CDS in dictonary speichern (PARALELLISIEREN)
-In short.tsv Gleichung für later_complete_ORFS berechnen und zu HC dictonary hinzufügen 
-In normal.tsv für incomplete_dictonary prüfen, ob C-terminus abgedeckt ist und Gleichung erfüllt, dann zu HC dictonary hinzufügen 
-
-Neue Transdecoder.pep file mit neuen Ergebnissen 
ODER alle complete durchgehen in bestehender pep file und bei incomplete
    schauen ob gekürzt wurde. Wenn ja zum nächsten Element.
'''

def load_config(config_file):
    with open(config_file, "r") as config_file:
        input_files = yaml.safe_load(config_file)
        return input_files

#MAIN
parser = argparse.ArgumentParser(description='Genome annotation with transcriptomic data like RNA-seq and Iso-seq data')  
parser.add_argument('-t', '--threads', default=4, help='Number of threads (default=4)', required=False)
parser.add_argument('-y', '--config', help='Config file input', metavar='<config.yaml>', required=False) #required=True

parser.add_argument('--isoseq', action='store_true', help='Use this option if you want to process isoseq data only')
parser.add_argument('--rnaseq', action='store_true', help='Use this option if you want to process rnaseq data only')
parser.add_argument('--mixed', action='store_true', help='Use this option if you want to process both rnaseq and isoseq data')

args = parser.parse_args()
threads = args.threads
input_files = load_config(args.config)
genome_file = input_files["genome"]
rnaseq_paired_sets = input_files.get("rnaseq_paired_sets", []) #Wenn Liste nicht vorhanden, dann leere Liste
rnaseq_single_sets = input_files.get("rnaseq_single_sets", [])
isoseq_sets = input_files.get("isoseq_sets", [])
protein_file = input_files["protein"] #optional?

#Intercept if given data doesnt match the chosen option  
if rnaseq_paired_sets == [] and rnaseq_single_sets == [] and isoseq_sets == []:
    print("Error: No transcriptomic data found in config file. Please provide at least one set of RNA-seq or Iso-seq data.")
    sys.exit(1)
if not args.rnaseq and not args.isoseq and not args.mixed:
    if rnaseq_paired_sets != [] or rnaseq_single_sets != []:
        args.rnaseq = True
    if isoseq_sets != []:
        args.isoseq = True
    print("You did not specify which data you want to process. The mode is set based on given data.")
if rnaseq_paired_sets == [] and rnaseq_single_sets == [] and args.rnaseq:
    print("Error: No RNA-seq data found in config file. Please provide at least one set of RNA-seq data.")
    sys.exit(1)
if isoseq_sets == [] and args.isoseq:
    print("Error: No Iso-seq data found in config file. Please provide at least one set of Iso-seq data.")
    sys.exit(1)
if (rnaseq_paired_sets == [] and rnaseq_single_sets == []) or (isoseq_sets == []) and args.mixed:
    print("Error: You chose the mixed option. Please provide both RNA-seq and Iso-seq data.")
    sys.exit(1)

process_rnaseq = args.rnaseq or args.mixed
process_isoseq = args.isoseq or args.mixed

'''
if process_rnaseq:
    #indexing(genome_file)
    alignments_list = mapping_short(rnaseq_paired_sets, rnaseq_single_sets)
    sam_to_bam(alignments_list)    
    if len(alignments_list) > 1:    #NOCHMAL PRÜFEN OB HIER WIRKLICH BAMFILES ÜBERGEBEN WERDEN
        alignment_rnaseq = file_name(merge_bam_files(alignments_list[0], alignments_list[1])) + ".bam"
    else:
        alignment_rnaseq = file_name(alignments_list[0]) + ".bam"

if process_isoseq:
    alignment_isoseq = mapping_long(genome_file, isoseq_sets)
    sam_file_list = [alignment_isoseq]
    sam_to_bam(sam_file_list) 
    alignment_isoseq = file_name(alignment_isoseq) + ".bam"
'''
#assembling(alignment_rnaseq, alignment_isoseq)  #Für alleine testen leer machen
#orfsearching(genome_file, "transcripts_mixed_test1.gtf")  #Vielleicht eher Was returned wurde als input übergeben
#protein_aligning(genome_file, protein_file, "/home/s-amknut/GALBA/tools/blosum62_1.csv") 
short_start_dict = shorten_incomplete_Orfs("transcripts_test1.fasta.transdecoder.pep")
#make_diamond_db(protein_file)
#validating_ORFs("shortened_candidates.pep", "diamond_shortened.tsv")
#make_diamond_db(protein_file)
#validating_ORFs("transcripts_test1.fasta.transdecoder.pep", "diamond_normal.tsv")
classifications_dict = get_cds_classification("diamond_normal.tsv", "diamond_shortened.tsv", short_start_dict)
classifications_hc_dict = get_optimized_pep_file("transcripts_test1.fasta.transdecoder.pep", "shortened_candidates.pep", classifications_dict, short_start_dict)

    

#make_diamond_db(protein_file)
#validating_ORFs("revised_candidates.pep", "diamond_revised.tsv")
#get_hc_cds("diamond_revised.tsv", classifications_hc_dict)

#transdecoder_id_dict = parse_transdecoder_file("revised_candidates.pep")
#from_transcript_to_genome_coords("transcripts_mixed_test1.gtf", transdecoder_id_dict)

#TO DOs:
#-Variablen und Funktionsnamen anpassen
#-Funktion die prüft ob files vorhanden wie CreateThis() von GeneMark 
#-cwd integrieren
#-verstehen was das --dta in hisat2 bedeutet 
#-prints überarbeiten
#-Ausgaben für subprocess.run überarbeiten, bei error mehr ausprinten lassen.
#-yaml-file prüfen, ob namen keine weiteren dots oder unterstriche enthalten sind. Alle Rnaseq files und alle isoseq files
# müssen dasselbe Format haben (Also unter sich)
#-Exception as e: print(e) einbauen
#-Vielleicht noch Option einbauen, dass nur einmal Pfad angegeben werden muss und sonst nur Namen der Files
#-Abfangen, wenn keine files in config liegen/Unter dem falschen Listennamen
#-Threads überall hinzufügen
#-f-strings da einfügen wo möglich
#-überlegen, wo Scoring Matrix eingefügt werden soll (Eine vorgeben oder von Nutzer hinzufügen lassen?)
#-Threads max. rausfinden und festlegen
#-Outputnamen vom Nutzer festlegen lassen und dann einzelne outputs an den Namen anpassen
#-Nutze seqIO Funktion um pep file inhalt überall in dictionary zu speichern: shortened_pep_dict = SeqIO.to_dict(SeqIO.parse(shortened_pep, "fasta"))
#-Optionen für einzelne Programe einführen, wie --genemark_path=/home/... oder --diamond_path=/home/...

#FRAGEN:
#-Sollte ich mit -G die stringtie Option nutzen, eine Referenzannotation zu verwenden? --> Diese dann in die yaml file oder parser?
#-Optionen in Ordnung oder noch Unterscheidung zwischen single und paired-end?
#-Ist es richtig, dass man single-end und paired-end beide nutzt? Oder wird in der Regel nur eins davon genutzt?
#-Soll die Übergabe von Proteinfiles optional sein?
#-Richtig, dass alle erstellten files in cwd gespeichert werden? Soll ich Funktion einfügen, dass man sich das aussuchen kann wohin?
#-Scoring Matrix von Nutzer einfügen lassen oder selbst eine vorgeben?
#-Ist es richtig, dass ein Trankript nach dem splicing mehrere CDS haben kann? D.h. dass man UTR mittig hat?

#-Wann Incomplete CDS als HC bezeichnet?
#-Training von AUGUSTUS auch in meiner Hand? Wenn ja, passiert das vorher?
#-Wieso steht in Doktorarbeit section 5.3.2.3 eine zweite Diamondsuche drin. Wir haben ja schon eine gemacht. 
#-Wie vergleiche ich die Annotationen. Welche ist die finale, die von Transdecoder?

#Plan:
#-Transdecoder macht ORF prediction -> Davor intron hints von spliced rnaseq daten wie bei genemark?
#-miniprot macht protein alignment
#-Proteinalignment von miniprot mit dem von Transdecoder vergleichen -> Prediction ergänzen oder verwerfen???
#-Mit miniprot trainingsgenen Augustus trainieren und hc hints an Augustus übergeben

#Oder:
#-Transdecoder macht ORF prediction -> Davor intron hints von spliced rnaseq daten wie bei genemark?
#-Diamond nutzt .pep file von Transdecoder und sucht homologe Proteine 
#-Spaln aligniert die homologen Proteine zurück ans Genom, um predictions genauer zu machen 

#In GeneMark:
#-3 Arten hints:
#-Transcript + protein support (Bei mir vielleicht Diamond)
#-Transcript + ab initio support (Bei mir vielleicht Augustus)
#-Nur Protein support (Bei mir vielleicht miniprot)

'''
Auswertung:
genemark.gtf mit allen rnaseq und protein.fa Daten:
#-----------------| Sensitivity | Precision  |  
        Base level:    94.1     |    91.2    |
        Exon level:    84.6     |    88.5    |
      Intron level:    90.5     |    90.3    |
Intron chain level:    55.4     |    75.6    |
  Transcript level:    57.7     |    77.1    |
       Locus level:    80.9     |    79.7    |

main.py mit allen rnaseq und isoseq Daten, nur bis transdecoder Aufruf (ohne Kategorisierung und ohne hc Filter):
#-----------------| Sensitivity | Precision  |
        Base level:    83.0     |    67.9    |
        Exon level:    60.9     |    57.3    |
      Intron level:    86.9     |    84.3    |
Intron chain level:    47.2     |    41.1    |
  Transcript level:    42.0     |    40.0    |
       Locus level:    55.5     |    65.9    |

main.py mit allen rnaseq und isoseq Daten, nur bis transdecoder Aufruf (ohne Kategorisierung und ohne hc Filter):
AUCH CDS VORHERGESAGT WENN NUR EIN EXON PRO TRANSCRIPT VORHANDEN

#-----------------| Sensitivity | Precision  |
        Base level:    83.2     |    65.3    |
        Exon level:    60.7     |    56.8    |
      Intron level:    86.8     |    82.5    |
Intron chain level:    42.5     |    36.6    |
  Transcript level:    37.5     |    35.7    |
       Locus level:    49.6     |    63.1    |

main.py mit allen rnaseq und isoseq Daten, nur bis transdecoder Aufruf (mit Kategorisierung und ohne hc Filter):

'''