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
    
#def correct_incomplete_Orfs(transdecoder_pep):
 #   try:
  #      with open("compare_cds.fasta", "w") as output:
   #         for record in SeqIO.parse(transdecoder_pep, "fasta"):
    #            if "type:5prime_partial" in record.description or "type:internal" in record.description:
     #               m_position = record.seq.find("M")
      #              #print("First Methionine in record ", record.id, " is located at Position ", m_position)
       #             if m_position == -1:
        #                record.id = record.id + ":none"
         #               SeqIO.write(record, output, "fasta")
          #          else:
           #             id = record.id
            #            record.id = id + ":long"
             #           SeqIO.write(record, output, "fasta")
              #          record.seq = record.seq[m_position:]
               #         record.id = id + ":short"
                #        SeqIO.write(record, output, "fasta")
               # else:
                #    record.id = record.id + ":none"
                 #   SeqIO.write(record, output, "fasta")
    
 #   except Exception:
  #      print("Could not write input file for Diamond.")
   #     sys.exit(1)

def correct_incomplete_Orfs(transdecoder_pep):
    try:
        short_start_dict = {}
        with open("shortened_candidates.pep", "w") as output:
            for record in SeqIO.parse(transdecoder_pep, "fasta"):
                if "type:5prime_partial" in record.description or "type:internal" in record.description:
                    m_position = record.seq.find("M")
                    if m_position != -1:
                        record.seq = record.seq[m_position:]
                        record.description = "Shortend sequence on from Position "+str(m_position) #NEU noch zu TESTEN (19.09.)
                        SeqIO.write(record, output, "fasta")
                        short_start_dict[record.id] = m_position
        return short_start_dict
    except Exception:
        print("Could not write input file for Diamond.")
        sys.exit(1)

def validating_ORFs(protein_file, shortened_pep, normal_pep):
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
    
    try:
        command = [
            "diamond",
            "blastp",
            "-d",
            "protein_db",
            "-q",
            shortened_pep,
            "-o",
            "diamond_shortened.tsv"
        ]
        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Blastp search completed successfully for complete CDS candidates")
        else:
            print("Error during blastp search for complete CDS candidates")
            print(result.stderr)
    
    except Exception:
        print("Could not run diamond blastp command for complete CDS candidates.")
        sys.exit(1)

    try:
        command = [
            "diamond",
            "blastp",
            "-d",
            "protein_db",
            "-q",
            normal_pep,
            "-o",
            "diamond_normal.tsv"
        ]
       # result = subprocess.run(command, capture_output=True)

        #if result.returncode == 0:
         #   print("Blastp search completed successfully for incomplete CDS candidates")
        #else:
         #   print("Error during blastp search for incomplete CDS candidates")
          #  print(result.stderr)
    
    except Exception:
        print("Could not run diamond blastp command for incomplete CDS candidates.")
        sys.exit(1)

def get_cds_classification(shortened_tsv, normal_tsv, short_start_dict):
    header_list = ["cdsID", "proteinID", "percIdentMatches", "alignLength", "mismatches", "gapOpenings", "alignStart", "alignEnd", "proteinStart", "proteinEnd", "eValue", "bitScore"]
    df_shortened = pd.read_csv(shortened_tsv, delimiter='\t', header=None, names=header_list)
    df_normal = pd.read_csv(normal_tsv, delimiter='\t', header=None, names=header_list)
    merged_df = pd.merge(df_shortened, df_normal, on=["cdsID", "proteinID"], suffixes=('_short', '_normal'))
    merged_df = merged_df.drop(columns=["alignLength_short", "alignLength_normal", "mismatches_short", "mismatches_normal", "gapOpenings_short", "gapOpenings_normal", "alignEnd_short", "alignEnd_normal", "proteinEnd_short", "proteinEnd_normal", "eValue_short", "eValue_normal"])
    merged_df["shortStart"] = merged_df["cdsID"].map(short_start_dict)
    #Muss noch dran denken die cdsIDs von denen es keine shorts gibt zu berücksichtigen, die sind nicht in merged_df
    #merged_df = merged_df.drop_duplicates(subset=["cdsID", "proteinID_short", "proteinID_normal"])
    merged_df["supportScore"] = None
    merged_df["classification"] = None
    #test_df = merged_df[merged_df.index<61]
    count = 0
    for i, cds in merged_df.iterrows():
    #for i, cds in test_df.iterrows():
        q_incomplete_start = cds["proteinStart_normal"]
        t_incomplete_start = cds["alignStart_normal"]
        t_complete_start = cds["alignStart_short"] + cds["shortStart"]  #not -1 because alignmentstart is 1-based but Mposition not
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
        if t_complete_start < t_incomplete_start:
            print("cdsID: ", cds["cdsID"] , "incomplete Start: ", t_incomplete_start, "complete Start: ", t_complete_start)
    
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
            print("t_complete_start: ", cds["alignStart_normal"], "t_incomplete_start: ", cds["alignStart_short"] )

    pd.set_option('display.max_columns', None)
    #print(merged_df.head())
    return classifications

def parse_transdecoder_file(transdecoder_pep):
    transdecoder_id_dict = {}
    for record in SeqIO.parse(transdecoder_pep, "fasta"):
        id = record.id
        id = id.split(".p")[0]
        #if id not in transdecoder_id_dict:
        #transdecoder_id_dict[id] = record.description.split(" ")
        describtion = record.description.split(" ")
        cds_transcript_coords = re.search(r":(\d+)-(\d+)\(\+\)", describtion[7])
        if cds_transcript_coords:
            start_cds_transcript = cds_transcript_coords.group(1)
            stop_cds_transcript = cds_transcript_coords.group(2)
            cds_length = int(stop_cds_transcript) - int(start_cds_transcript) + 1
            pair = (start_cds_transcript, stop_cds_transcript, cds_length)
            if id not in transdecoder_id_dict:
                transdecoder_id_dict[id] = [pair]
            else:
                transdecoder_id_dict[id].append(pair)

    #print(transdecoder_id_dict)
    return transdecoder_id_dict

'''
        more_cds = True
        cds_number = 1
        while more_cds:
            transcript_id_transdecoder = transcript_id + ".p" + str(cds_number)
            if transcript_id_transdecoder in transdecoder_id_dict:
                describtion = transdecoder_id_dict[transcript_id_transdecoder]
                cds_transcript_coords = re.search(r":(\d+)-(\d+)\(\+\)", describtion[7])
                if cds_transcript_coords:
                    start_cds_transcript = cds_transcript_coords.group(1)
                    stop_cds_transcript = cds_transcript_coords.group(2)
                    print("ID: ", transcript_id_transdecoder, "Start in Transkript: ", start_cds_transcript, "Stop in Transkript: ", stop_cds_transcript)
                    start_cds_genome = int(start_genome) + int(start_cds_transcript) - 1
                    stop_cds_genome = int(start_genome) + int(stop_cds_transcript) - 1
                    output.write(f"{seqname}\tPreGalba\tcds\t{start_cds_genome}\t{stop_cds_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                cds_number += 1
        #print(transdecoder_id_dict)
        return transdecoder_id_dict
'''
 
def in_both_dicts(dict1, dict2):
    not_in_dict1 = {}
    count = 0
    for key in dict1:
        if key not in dict2:
            not_in_dict1[key] = dict1[key]
            count +=1
            #print(key)
    #print(not_in_dict1)
    #print("Id in Trans but not in Stringtie: ", count)

def from_transcript_to_genome_coords(stringtie_gtf, transdecoder_pep, transdecoder_id_dict):
    annotation_file = "annotation.gtf"
    test = 0
    with open(transdecoder_pep, "r") as transdecoder, open(stringtie_gtf, "r") as stringtie, open(annotation_file, "w") as output:
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
                    #if test == 1:
                     #   break
                    #test = 1
                    '''
                    print("Listenlänge: ", len(exon_coords_list))
                    if len(exon_coords_list) > 0:
                        cds_index = 0
                        exon_index = 0
                        last_transcript_id = exon_coords_list[0][2]
                        #Macht keinen sinn weil schon die transcript id vom nächsten Transkript genommen wird
                        if last_transcript_id in transdecoder_id_dict:
                            cds_list = transdecoder_id_dict[last_transcript_id]
                            cds_start_transcript = cds_list[cds_index][0]
                            cds_stop_transcript = cds_list[cds_index][1]
                            cds_start_genome = int(start_genome) + int(cds_start_transcript) - 1
                            cds_stop_genome = int(start_genome) + int(cds_stop_transcript) - 1
                            five_prime_UTR_stop_transcript = int(cds_start_transcript) - 1
                            five_prime_UTR_stop_genome = int(start_genome) + five_prime_UTR_stop_transcript - 1
                            output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{start_genome}\t{five_prime_UTR_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{last_transcript_id}\";\n")
                            for exon in exon_coords_list:
                                exon_start_genome = exon[0]
                                exon_stop_genome = exon[1]
                                if cds_stop_genome > int(exon_stop_genome):
                                    cds_stop_genome = int(exon_stop_genome)
                                    output.write(f"{seqname}\tPreGalba\tcds\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                    cds_start_genome = cds_stop_genome + 1
                                else:
                                    output.write(f"{seqname}\tPreGalba\tcds\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                        exon_coords_list.clear()
                    '''
                        
                    output.write(f"{seqname}\tPreGalba\tgenome\t{start_genome}\t{stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    output.write(f"{seqname}\tPreGalba\ttranscript\t{start_genome}\t{stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    if transcript_id in transdecoder_id_dict:
                        cds_index = 0
                        cds_current_length = 0
                        new_cds = True
                    else:
                        new_cds = False
                    '''
                    more_cds = True
                    cds_number = 1
                    while more_cds:
                        transcript_id_transdecoder = transcript_id + ".p" + str(cds_number)
                        if transcript_id_transdecoder in transdecoder_id_dict:
                            describtion = transdecoder_id_dict[transcript_id_transdecoder]
                            cds_transcript_coords = re.search(r":(\d+)-(\d+)\(\+\)", describtion[7])
                            if cds_transcript_coords:
                                start_cds_transcript = cds_transcript_coords.group(1)
                                stop_cds_transcript = cds_transcript_coords.group(2)
                                print("ID: ", transcript_id_transdecoder, "Start in Transkript: ", start_cds_transcript, "Stop in Transkript: ", stop_cds_transcript)
                                start_cds_genome = int(start_genome) + int(start_cds_transcript) - 1
                                stop_cds_genome = int(start_genome) + int(stop_cds_transcript) - 1
                                output.write(f"{seqname}\tPreGalba\tcds\t{start_cds_genome}\t{stop_cds_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                            cds_number += 1
                        #else:
                            #5'UTR und 3'UTR noch hinzufügen
                         #   more_cds = False 
                '''
                if feature == 'exon':
                    output.write(f"{seqname}\tPreGalba\texon\t{start_genome}\t{stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    exon_number = re.search(r'exon_number "([^"]+)"', attributes)
                    exon_number = exon_number.group(1)
                    #exon_coords_list.append((start_genome, stop_genome, transcript_id))
                    
                    if exon_number == "1":
                        previous_exon_stop = stop_genome
                       # exon_coords_list.append((start_genome, stop_genome, transcript_id))
                        #new_cds = True
                        #cds_index = 0
                    else:
                        intron_start_genome = int(previous_exon_stop) + 1
                        intron_stop_genome = int(start_genome) - 1
                        previous_exon_stop = stop_genome
                        output.write(f"{seqname}\tPreGalba\tintron\t{intron_start_genome}\t{intron_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
        
                            
                    if transcript_id in transdecoder_id_dict:
                        cds_list = transdecoder_id_dict[transcript_id]
                        if new_cds == True:
                            if len(cds_list) > cds_index:
                                cds_start_transcript = cds_list[cds_index][0]
                                cds_stop_transcript = cds_list[cds_index][1]
                                cds_total_length = cds_list[cds_index][2]
                                print("cds total length: ", cds_total_length)
                                cds_start_genome = int(start_genome) + int(cds_start_transcript) - 1
                                cds_stop_genome = int(start_genome) + int(cds_stop_transcript) - 1
                                new_cds = False
                            else:
                                continue
                       # else:
                        #    cds_start_genome = cds_overstop + intron_stop_genome - intron_start_genome + 2
                        print("cds_cuurent_length: ", cds_current_length)
                        #if cds_stop_genome > int(stop_genome):
                        if cds_current_length+int(stop_genome)-int(start_genome)+1<cds_total_length:
                            print("hi")
                            if cds_index != 0:
                                print("CDS Index: ", cds_index)
                                cds_start_genome = int(start_genome)
                            output.write(f"{seqname}\tPreGalba\tcds\t{cds_start_genome}\t{stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                            cds_current_length += int(stop_genome) - cds_start_genome + 1
                            #cds_overstop = int(stop_genome)
                        else:
                           # output.write(f"{seqname}\tPreGalba\tcds\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                            #new_cds = True
                            #cds_index += 1
                            output.write(f"{seqname}\tPreGalba\tcds\t{start_genome}\t{cds_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                            new_cds = True
                            cds_index += 1
                        

                        #Idee: In der Mitte einfach exon = cds und nur für erstes und letztes exon UTRs bestimmen

                    '''
                    if cds_index == 0:
                        cds_list = transdecoder_id_dict[transcript_id]
                        cds_start_transcript = cds_list[cds_index][0]
                        cds_stop_transcript = cds_list[cds_index][1]
                        five_prime_UTR_stop_transcript = int(cds_start_transcript) - 1
                        five_prime_UTR_stop_genome = int(start_genome) + five_prime_UTR_stop_transcript - 1
                        output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{start_genome}\t{five_prime_UTR_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                        cds_start_genome = int(start_genome) + int(cds_start_transcript) - 1
                        cds_stop_genome = int(start_genome) + int(cds_stop_transcript) - 1
                        if cds_stop_genome > int(stop_genome):
                            cds_stop_genome = int(stop_genome)
                            output.write(f"{seqname}\tPreGalba\tcds\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                        else:
                            output.write(f"{seqname}\tPreGalba\tcds\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    '''

                        #start_cds_genome = int(start_genome) + int(start_cds_transcript) - 1
                        #stop_cds_genome = int(start_genome) + int(stop_cds_transcript) - 1
                        #output.write(f"{seqname}\tPreGalba\tcds\t{start_cds_genome}\t{stop_cds_genome}\t.\t{strand}\t{frame}\tgene_id \"{gene_id}\"
                        #cds_list = transdecoder_id_dict[transcript_id]
                        #if cds_list[cds_index][0] 
                        #print(cds_list)
                    



        #in_both_dicts(transdecoder_id_dict, stringtie_dict)         
'''
Chr1	StringTie	transcript	3676	5861	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "27.170204"; FPKM "4.147030"; TPM "5.471606";
Chr1	StringTie	exon	3676	3913	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "19.891108";

>STRG.1.1.p1 GENE.STRG.1.1~~STRG.1.1.p1  ORF type:5prime_partial (+),score=23.30 len:152 STRG.1.1:3-461(+)
QSRQRNSGSYNTYSEYDSANHGQQFNENSNIMQQQPLQGSFNPLLEYDFANHGGQWLSDY

Chr1    MAKER   transcript      1000    4500    .       +       .       gene_id "gene1"; transcript_id "transcript2"; gene_name "Gene1"; 
Chr1    MAKER   exon            1000    1500    .       +       .       gene_id "gene1"; transcript_id "transcript2"; exon_number "1"; 
Chr1    MAKER   exon            2500    3500    .       +       .       gene_id "gene1"; transcript_id "transcript2"; exon_number "2"; 
Chr1    MAKER   exon            4000    4500    .       +       .       gene_id "gene1"; transcript_id "transcript2"; exon_number "3"; 
Chr1    MAKER   CDS             1100    1500    .       +       0       gene_id "gene1"; transcript_id "transcript2"; 
Chr1    MAKER   CDS             2500    3500    .       +       2       gene_id "gene1"; transcript_id "transcript2"; 
Chr1    MAKER   CDS             4000    4500    .       +       0       gene_id "gene1"; transcript_id "transcript2"; 
Chr1    MAKER   five_prime_UTR  1000    1099    .       +       .       gene_id "gene1"; transcript_id "transcript2"; 
Chr1    MAKER   three_prime_UTR 4501    4500    .       +       .       gene_id "gene1"; transcript_id "transcript2"; 

Chr1    Araport11       CDS     3760    3913    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; cds_type "Initial"; count "1_6";
Chr1    Araport11       CDS     3996    4276    .       +       2       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; cds_type "Internal"; count "2_6";
Chr1    Araport11       CDS     4486    4605    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; cds_type "Internal"; count "3_6";
Chr1    Araport11       CDS     4706    5095    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; cds_type "Internal"; count "4_6";
Chr1    Araport11       CDS     5174    5326    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; cds_type "Internal"; count "5_6";
Chr1    Araport11       CDS     5439    5630    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; cds_type "Terminal"; count "6_6";
Chr1    Araport11       intron  3914    3995    .       +       1       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; count "1_5"; site_seq "GT_AG";
Chr1    Araport11       intron  4277    4485    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; count "2_5"; site_seq "GT_AG";
Chr1    Araport11       intron  4606    4705    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; count "3_5"; site_seq "GT_AG";
Chr1    Araport11       intron  5096    5173    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; count "4_5"; site_seq "GT_AG";
'''

#-Wann Incomplete CDS als HC bezeichnet?
#-Wie kann ich aus Kategorisierung, also aus file mit neuen cds sequenzen, die GFF3/GTF Datei erstellen?

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
            > GMS-T log-odds score > 50,
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
'''
Normal:
STRG.10016.1.p1	prot521	56.0	125	55	0	1	125	187	311	7.82e-43	143
STRG.10016.1.p1	prot853	35.8	123	70	4	8	125	159	277	2.29e-15	67.0
STRG.10016.1.p1	prot852	36.8	117	65	5	8	119	143	255	5.01e-14	63.2
STRG.10017.1.p2	prot660	29.4	68	42	2	1	64	61	126	6.42e-05	38.1
STRG.10034.1.p1	prot1315	47.4	137	40	4	34	169	20	125	7.81e-26	103
STRG.10039.1.p1	prot1288	56.7	127	49	4	40	162	28	152	1.66e-40	139
STRG.10039.1.p1	prot837	35.6	132	79	2	37	162	21	152	8.57e-23	89.7
STRG.10057.1.p1	prot1266	33.8	195	59	9	65	256	43	170	3.02e-18	76.3
STRG.10090.1.p1	prot691	28.2	387	231	8	17	382	66	426	3.35e-39	143
STRG.10101.1.p1	prot64	56.2	169	69	3	1	168	146	310	1.65e-54	175
STRG.10101.1.p3	prot64	58.1	31	13	0	47	77	89	119	1.19e-04	35.8 

Shortened:
STRG.10016.1.p1	prot521	55.8	113	50	0	2	114	199	311	1.79e-37	128
STRG.10016.1.p1	prot853	35.6	118	67	4	2	114	164	277	3.55e-14	63.2
STRG.10016.1.p1	prot852	36.6	112	62	5	2	108	148	255	5.75e-13	59.7
STRG.10039.1.p1	prot1288 56.7	127	49	4	29	151	28	152	1.14e-40	139
STRG.10039.1.p1	prot837	35.6	132	79	2	26	151	21	152	6.26e-23	89.7
STRG.10110.1.p1	prot1201 44.4	99	55	0	148	246	803	901	1.01e-22	92.8
STRG.10174.1.p1	prot1095 40.7	216	124	3	28	242	316	528	7.82e-50	169
STRG.10174.1.p1	prot446	31.4	242	162	3	1	242	609	846	9.54e-45	156
STRG.10174.1.p1	prot447	31.4	242	162	3	1	242	609	846	9.54e-45	156
STRG.10174.1.p1	prot448	31.4	242	162	3	1	242	609	846	9.54e-45	156

STRINGTIE:
Chr1	StringTie	transcript	3676	5861	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "27.170204"; FPKM "4.147030"; TPM "5.471606";
Chr1	StringTie	exon	3676	3913	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "19.891108";
Chr1	StringTie	exon	3996	4276	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "15.570964";
Chr1	StringTie	exon	4486	4605	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "3"; cov "22.015558";
Chr1	StringTie	exon	4706	5095	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "4"; cov "31.399448";
Chr1	StringTie	exon	5174	5326	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "5"; cov "34.461781";
Chr1	StringTie	exon	5439	5861	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "6"; cov "33.396362";

STRG.1.1	transdecoder	gene	1	1157	.	+	.	ID=GENE.STRG.1.1~~STRG.1.1.p1;Name="ORF type:5prime_partial (+),score=23.30"
STRG.1.1	transdecoder	mRNA	1	1157	.	+	.	ID=STRG.1.1.p1;Parent=GENE.STRG.1.1~~STRG.1.1.p1;Name="ORF type:5prime_partial (+),score=23.30"
STRG.1.1	transdecoder	exon	1	1157	.	+	.	ID=STRG.1.1.p1.exon1;Parent=STRG.1.1.p1
STRG.1.1	transdecoder	CDS   	3	461	    .	+	0	ID=cds.STRG.1.1.p1;Parent=STRG.1.1.p1
STRG.1.1	transdecoder	three_prime_UTR	462	1157	.	+	.	ID=STRG.1.1.p1.utr3p1;Parent=STRG.1.1.p1

>STRG.1.1.p1 GENE.STRG.1.1~~STRG.1.1.p1  ORF type:5prime_partial (+),score=23.30 len:152 STRG.1.1:3-461(+)
QSRQRNSGSYNTYSEYDSANHGQQFNENSNIMQQQPLQGSFNPLLEYDFANHGGQWLSDY

Chr1	StringTie	transcript	6812	8720	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; cov "171.187973"; FPKM "5.001022"; TPM "7.711779";
Chr1	StringTie	exon	6812	7069	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; cov "125.622360";
Chr1	StringTie	exon	7157	7232	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "2"; cov "230.052353";
Chr1	StringTie	exon	7384	7450	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "3"; cov "211.434280";
Chr1	StringTie	exon	7564	7649	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "4"; cov "219.834091";
Chr1	StringTie	exon	7762	7835	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "5"; cov "211.385376";
Chr1	StringTie	exon	7942	7987	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "6"; cov "209.836029";
Chr1	StringTie	exon	8236	8325	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "7"; cov "222.475906";
Chr1	StringTie	exon	8417	8464	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "8"; cov "227.134079";
Chr1	StringTie	exon	8571	8720	1000	-	.	gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "9"; cov "93.510895";

STRG.9450.3	transdecoder	gene	1	1609	.	+	.	ID=GENE.STRG.9450.3~~STRG.9450.3.p1;Name="ORF type:complete (+),score=47.10"
STRG.9450.3	transdecoder	mRNA	1	1609	.	+	.	ID=STRG.9450.3.p1;Parent=GENE.STRG.9450.3~~STRG.9450.3.p1;Name="ORF type:complete (+),score=47.10"
STRG.9450.3	transdecoder	five_prime_UTR	1	840	.	+	.	ID=STRG.9450.3.p1.utr5p1;Parent=STRG.9450.3.p1
STRG.9450.3	transdecoder	exon	1	1609	.	+	.	ID=STRG.9450.3.p1.exon1;Parent=STRG.9450.3.p1
STRG.9450.3	transdecoder	CDS	841	1506	.	+	0	ID=cds.STRG.9450.3.p1;Parent=STRG.9450.3.p1
STRG.9450.3	transdecoder	three_prime_UTR	1507	1609	.	+	.	ID=STRG.9450.3.p1.utr3p1;Parent=STRG.9450.3.p1

ZIEL:
Chr1    MAKER   gene            1000    5000    .       +       .       gene_id "gene1"; gene_name "Gene1"; 
Chr1    MAKER   transcript      1000    5000    .       +       .       gene_id "gene1"; transcript_id "transcript1"; gene_name "Gene1"; 

# Exons for the first transcript
Chr1    MAKER   exon            1000    1200    .       +       .       gene_id "gene1"; transcript_id "transcript1"; exon_number "1"; 
Chr1    MAKER   exon            2000    3000    .       +       .       gene_id "gene1"; transcript_id "transcript1"; exon_number "2"; 
Chr1    MAKER   exon            4000    5000    .       +       .       gene_id "gene1"; transcript_id "transcript1"; exon_number "3"; 

# Coding sequence (CDS) within these exons
Chr1    MAKER   CDS             1050    1200    .       +       0       gene_id "gene1"; transcript_id "transcript1"; 
Chr1    MAKER   CDS             2000    3000    .       +       2       gene_id "gene1"; transcript_id "transcript1"; 
Chr1    MAKER   CDS             4000    4500    .       +       1       gene_id "gene1"; transcript_id "transcript1"; 

# UTRs for completeness
Chr1    MAKER   five_prime_UTR  1000    1049    .       +       .       gene_id "gene1"; transcript_id "transcript1"; 
Chr1    MAKER   three_prime_UTR 4501    5000    .       +       .       gene_id "gene1"; transcript_id "transcript1"; 

# A second transcript of the same gene with slightly different structure
Chr1    MAKER   transcript      1000    4500    .       +       .       gene_id "gene1"; transcript_id "transcript2"; gene_name "Gene1"; 
Chr1    MAKER   exon            1000    1500    .       +       .       gene_id "gene1"; transcript_id "transcript2"; exon_number "1"; 
Chr1    MAKER   exon            2500    3500    .       +       .       gene_id "gene1"; transcript_id "transcript2"; exon_number "2"; 
Chr1    MAKER   exon            4000    4500    .       +       .       gene_id "gene1"; transcript_id "transcript2"; exon_number "3"; 
Chr1    MAKER   CDS             1100    1500    .       +       0       gene_id "gene1"; transcript_id "transcript2"; 
Chr1    MAKER   CDS             2500    3500    .       +       2       gene_id "gene1"; transcript_id "transcript2"; 
Chr1    MAKER   CDS             4000    4500    .       +       0       gene_id "gene1"; transcript_id "transcript2"; 
Chr1    MAKER   five_prime_UTR  1000    1099    .       +       .       gene_id "gene1"; transcript_id "transcript2"; 
Chr1    MAKER   three_prime_UTR 4501    4500    .       +       .       gene_id "gene1"; transcript_id "transcript2"; 

'''
def compare_annotation(transdecoder_gff3, example_gff3):
    try:
        command = [
            "/home/s-amknut/GALBA/tools/gffcompare",
            "-r",
            example_gff3,
            "-o comparison",
            transdecoder_gff3
        ]
        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Comparison completed successfully")
        else:
            print("Error during comparison")
            print(result.stderr)
    
    except Exception:
        print("Could not run gffcompare command.")
        print(result.stderr)
        print(result.stdout)
        sys.exit(1)


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
#short_start_dict = correct_incomplete_Orfs("transcripts.fasta.transdecoder.pep")
#validating_ORFs(protein_file, "shortened_candidates.pep", "transcripts.fasta.transdecoder.pep")
#classifications_dict = get_cds_classification("diamond_shortened.tsv", "diamond_normal.tsv", short_start_dict)
#annot = "/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/annot/pseudo.gff3"
#compare_annotation("transdecoder.fasta.transdecoder.gff3", annot)
#transdecoder_id_dict = parse_transdecoder_file("transcripts.fasta.transdecoder.pep")
#from_transcript_to_genome_coords("transcripts.gtf", "transcripts.fasta.transdecoder.cds", transdecoder_id_dict)
#transdecoder_id_dict = parse_transdecoder_file("transcripts_test1.fasta.transdecoder.pep")
#from_transcript_to_genome_coords("transcripts_mixed_test1.gtf", "transcripts_test1.fasta.transdecoder.pep", transdecoder_id_dict)
transdecoder_id_dict = parse_transdecoder_file("transcripts.fasta.transdecoder.pep")
from_transcript_to_genome_coords("transcripts.gtf", "transcripts.fasta.transdecoder.pep", transdecoder_id_dict)




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