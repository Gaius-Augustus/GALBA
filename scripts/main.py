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
from Bio.SeqRecord import SeqRecord

#FUNCTIONS
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
    try:
        hisat2_build_command = [
            "hisat2-build", 
            "--quiet",
            "-p",
            str(threads),        
            genome_fasta,           
            "genome"   
        ]

        print("Building genome index...")
    
        result = subprocess.run(hisat2_build_command, capture_output=False)
    
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
            output_1 = "alignment_single_rnaseq.sam" 
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
            output_2 = "alignment_paired_rnaseq.sam"  
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
        output_sam = "alignment_isoseq.sam" #threads neu prüfen
        # vorher: minimap2_command = ["minimap2", "-ax", "splice", "-uf", "-C5", genome, "-t", str(threads)] + isoseq_sets + ["-o", output_sam]
        minimap2_command = ["minimap2", "-ax", "splice:hq", "-uf", genome, "-t", str(threads)] + isoseq_sets + ["-o", output_sam]
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
        output_bam = "alignment_merged_rnaseq.bam"
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
        output_gtf = "transcripts.gtf"
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
        
        output_fa = "transcripts.fasta" 
        '''
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
        '''
        
        fasta_file_command = [
            "/opt/TransDecoder/util/gtf_genome_to_cdna_fasta.pl",
            transcripts_gtf,
            genome_fa,
        ]
        print("Constructing the transcript fasta file using the genome and the transcripts.gtf")
        print("Running command:", "".join(fasta_file_command))
        with open("transcripts.fasta", "w") as output:
            result = subprocess.run(fasta_file_command, stdout=output, stderr=subprocess.PIPE)

        if result.returncode == 0:
            print("Transcripts prepared successfully")
        else:
            print("Error during preparing transcripts")
            print(result.stderr)    
    
    except Exception: #ÄNDERN
        print("Could not run gffread command.")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    #No threads option for Transdecoder (ChatGPT says option would be to split input files)
    try:
        print("Extract the long open reading frames...")
        longORF_command = [
            "TransDecoder.LongOrfs",
            #"-S",
            "-t",
            "transcripts.fasta"
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
            "transcripts.fasta"
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

def convert_gtf_to_gff3(transcripts_gtf):
    try:
        command = [
            "/opt/TransDecoder/util/gtf_to_alignment_gff3.pl",
            transcripts_gtf,
            #">",
            #"transcripts.gff3"
        ]
    
        with open("transcripts.gff3", "w") as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE)

        if result.returncode == 0:
            print("Converted GTF to GFF3 file completed successfully.")

        else:
            print("Error during converting GTF to GFF3 file.")
    
    except Exception:
        print("Could not run util/gtf_to_alignment_gff3.pl module.") 
'''
def orfsearching(genome_fa, transcripts_gtf):
    output_fa = "transcripts_util.fasta"
    
    try:
        output_fa = "transcripts.fasta" 
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
'''
def protein_aligning(genome, protein, alignment_scoring):
    try: 
        command = [
            "/home/s-amknut/GALBA/tools/miniprot/miniprot",
           # "miniprot",
            "-t",       #NEU:threads
            str(threads),
            "--genome",
            genome,
            "--protein",
            protein,
            "--aln",
        ]
        print("Aligning proteins to genome...")
        with open("miniprot.aln", "w") as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE)
        if result.returncode == 0:
            print("Proteins aligned successfully")
        else:
            print("Error during protein alignment with miniprot")
            print("stdout:", result.stdout.decode())
            print("stderr:", result.stderr.decode())

    except Exception:
        print("Could not run miniprot command.")
        print(result.stderr)
        sys.exit(1)

    try: 
        #NEU: geändert von fstring zu liste
        command = f"/home/s-amknut/GALBA/tools/miniprot-boundary-scorer/miniprot_boundary_scorer -o miniprot_parsed.gff -s {alignment_scoring} < miniprot.aln"
        #command = [
         #   "/home/s-amknut/GALBA/tools/miniprot-boundary-scorer/miniprot_boundary_scorer",
          #  "-o",
           # "miniprot_parsed.gff",
            #"-s",
            #alignment_scoring,
            #"<",
            #"miniprot.aln"
        #]
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
                        new_length = len(record.seq) - 1 #damit stopp nicht mitgezählt wird
                        if coords:
                            old_start = int(coords.group(1)) 
                            new_start = old_start + m_position*3 #Nicht -1, denn m_position ist 0-based
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
            print("Blastp search completed successfully and ", output_tsv, " was created.")
        else:
            print("Error during diamond blastp search.")
            print(result.stderr)
    
    except Exception:
        print("Could not run diamond blastp command.")
        sys.exit(1)

def get_cds_classification(normal_tsv, shortened_tsv, short_start_dict):
    #Weil hier nach CDS-ID und nach Protein gemerged wird, werden einzelne herausgefiltert, die nicht in beiden vorkommen. 
    header_list = ["cdsID", "proteinID", "percIdentMatches", "alignLength", "mismatches", "gapOpenings", "queryStart", "queryEnd", "targetStart", "targetEnd", "eValue", "bitScore"]
    df_shortened = pd.read_csv(shortened_tsv, delimiter='\t', header=None, names=header_list)
    df_normal = pd.read_csv(normal_tsv, delimiter='\t', header=None, names=header_list)
    merged_df = pd.merge(df_shortened, df_normal, on=["cdsID", "proteinID"], suffixes=('_short', '_normal'))
    merged_df = merged_df.drop(columns=["alignLength_short", "alignLength_normal", "mismatches_short", "mismatches_normal", "gapOpenings_short", "gapOpenings_normal", "queryEnd_short", "queryEnd_normal", "targetEnd_short", "targetEnd_normal", "eValue_short", "eValue_normal"])
    merged_df["shortStart"] = merged_df["cdsID"].map(short_start_dict)
    merged_df["supportScore"] = None
    #merged_df["classification"] = None
    count = 0 
    
    for i, cds in merged_df.iterrows():
        if count >= 0 :
            count += 1
            #print(cds["shortStart"])
            q_incomplete_start = cds["queryStart_normal"]
            t_incomplete_start = cds["targetStart_normal"]
            t_complete_start = cds["targetStart_short"] #+ cds["shortStart"]  doch nicht weil eh in Protein#not -1 because alignmentstart is 1-based but Mposition not. +shortstart, weil Differenz von normal zu short gebraucht wird
            aai_incomplete = cds["percIdentMatches_normal"] 
            aai_complete = cds["percIdentMatches_short"] 
            if aai_complete == 0:
                aai_complete = 0.0001
            #if aai_incomplete == 0:
            #   aai_incomplete = 0.0001
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
            classifications[cds["cdsID"]] = "incomplete" #112 13.8% NEU 354 
        else:
            classifications[cds["cdsID"]] = "complete" #699 86.2% NEU 4714 
        
    #18.11.:
    #Es gab in shortened.pep file 5192 Einträge 
    #4714 werden gekürzt 90.8%
    #478 bleiben incomplete 9.2%
    #24.11.: 5068 candidates
    #1758 bleiben incomplete 
    pd.set_option('display.max_columns', None)
    #print(merged_df.head())
    countc = 0
    for id in classifications:
        if classifications[id] == "incomplete":
            countc += 1

    print("Incomplete: ", countc)
    return classifications
   
#Noch die Klassifkation rausnehmen
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
                if id in classifications: #Falls es Protein gab, das aligniert hat
                    if classifications[id] == "incomplete":
                        SeqIO.write(record, output, "fasta")  
                    else:
                        seq = shortened_pep_dict[id][0]
                        record.seq = seq
                        description = shortened_pep_dict[id][1]
                        #description = re.sub(r"len:\d+", f"len:{len(seq)-1}", record.description) #-1 damit * nicht gezählt wird
                        #coords = re.search(r":(\d+)-(\d+)\([\+\-]\)", description)
                        #start_normal = int(coords.group(1))
                        #start_short = start_normal + short_start_dict[id] * 3 #Nicht -1, denn m_position ist 0-based
                        #stop = int(coords.group(2))
                        #strand = coords.group(0)[-2]
                        #new_coords = f'{start_short}-{stop}({strand})'
                        #description = re.sub(r'(\d+-\d+\([\+\-]\))', new_coords, description)
                        if "type:5prime_partial" in description:
                            description = description.replace("type:5prime_partial", "type:complete")
                            classification = "complete"
                        else:
                            description = description.replace("type:internal", "type:3prime_partial")
                            classification = "3prime_partial"
                        record.description = description   
                        SeqIO.write(record, output, "fasta")
                else:
                    SeqIO.write(record, output, "fasta") #If there is no classification, there was no protein evidence found for the CDS 
                    #classifications[id] = "incomplete"   #(drüber nachdenken ob es sein kann, dass es nur für den normalen aber dafür für den kurzen evidenz gibt. Auch dann kommt cds in merged_df nicht vor)
            else:
                SeqIO.write(record, output, "fasta")

            #if classification == "complete" or classification == "5prime_partial":
             #   classifications_for_hc[id] = [classification, record.description]

    #return classifications_for_hc

    #4691 sind in classifications nicht drin, aber in normal_pep als 5prime/ internal
    #5502 sind in 5prime/internal in normal_pep
    #310 sind in normal_pep als 5prime_partial oder internal, aber nicht in short_pep_dict, weil sie kein M enthalten!!!
    #4682 sind in normal_pep als 5prime_partial oder internal, aber nicht in short_tsv
    #820 elemente sind in short_tsv
    #9 Elemente sind in short_tsv aber nicht in classifications -> Vermutung: Protein aligniert nur mit short oder nur mit normal und alignment taucht damit nicht in mergeddf auf

def from_pep_file_to_gff3(orf_pep, transcript_gtf, output_name):
    transcript_lengths = {}
    transcript_length = 0
    with open(transcript_gtf, "r") as transcript_file:
        for line in transcript_file:
            if line.startswith("#"):
                continue
            else:
                part = line.strip().split('\t')
                if "exon" == part[2]:
                    start = int(part[3])
                    stop = int(part[4])
                    length = stop - start + 1 
                    transcript_length += length
                    transcript_id = re.search(r'transcript_id "([^"]+)"', part[8]).group(1)
                    transcript_lengths[transcript_id] = transcript_length
                else:
                    transcript_length = 0
    with open(output_name, "w") as output:
        for record in SeqIO.parse(orf_pep, "fasta"):
            id_transdecoder = record.id
            id_stringtie = id_transdecoder.split(".p")
            id_stringtie = id_stringtie[0]
            tool = "transdecoder"
            description = record.description
            coords = re.search(r":(\d+)-(\d+)\([\+\-]\)", description)
            orf_start = int(coords.group(1))
            orf_stop = int(coords.group(2))
            transcript_length = transcript_lengths[id_stringtie]
            strand = re.search(r"\((\+|-)\)", description) 
            strand = strand.group(1) 
            description_parts = description.split()
            gene_id = description_parts[1]
            gene_name = description_parts[2] + " "+ description_parts[3] +" "+ description_parts[4]  
            output.write(f"{id_stringtie}\t{tool}\tgene\t1\t{transcript_length}\t.\t{strand}\t.\tID={gene_id};Name=\"{gene_name}\"\n")
            output.write(f"{id_stringtie}\t{tool}\tmRNA\t1\t{transcript_length}\t.\t{strand}\t.\tID={id_transdecoder};Parent={gene_id};Name=\"{gene_name}\"\n")
            if strand == "+" and orf_start > 1:
                output.write(f"{id_stringtie}\t{tool}\tfive_prime_UTR\t1\t{orf_start-1}\t.\t{strand}\t.\tID={id_transdecoder}.utr5p1;Parent={id_transdecoder}\n")
            if strand == "-" and orf_stop < transcript_length:
                output.write(f"{id_stringtie}\t{tool}\tfive_prime_UTR\t{orf_stop+1}\t{transcript_length}\t.\t{strand}\t.\tID={id_transdecoder}.utr5p1;Parent={id_transdecoder}\n")
            output.write(f"{id_stringtie}\t{tool}\texon\t1\t{transcript_length}\t.\t{strand}\t.\tID={id_transdecoder}.exon1;Parent={id_transdecoder}\n")
            output.write(f"{id_stringtie}\t{tool}\tCDS\t{orf_start}\t{orf_stop}\t.\t{strand}\t0\tID=cds.{id_transdecoder};Parent={id_transdecoder}\n")
            if strand == "-" and orf_start > 1:
                output.write(f"{id_stringtie}\t{tool}\tthree_prime_UTR\t1\t{orf_start-1}\t.\t{strand}\t.\tID={id_transdecoder}.utr3p1;Parent={id_transdecoder}\n")
            if strand == "+" and orf_stop < transcript_length:
                output.write(f"{id_stringtie}\t{tool}\tthree_prime_UTR\t{orf_stop+1}\t{transcript_length}\t.\t{strand}\t.\tID={id_transdecoder}.utr3p1;Parent={id_transdecoder}\n")
            output.write("\n")

def preparing_miniprot_gff_for_conflict_comparison(miniprot_gff):
    with open(miniprot_gff, "r") as gff_file, open("reference.bed", "w") as output:
        for line in gff_file: 
            if "CDS" in line:
                parts = line.split("\t")
                chromosome = parts[0]
                start = parts[3]
                stop = parts[4]
                prot_id = re.search(r"prot=(prot\d+)", line).group(1)
                strand = parts[6]
                output.write(f"{chromosome}\t{start}\t{stop}\t{prot_id}\t.\t{strand}\n")

def from_dict_to_pep_file(input_dict, output_name):
    with open(output_name, "w") as output:
        for entry in input_dict:
            # Create a SeqRecord directly using the Seq and description from the dictionary
            record = SeqRecord(
                input_dict[entry][1],  
                id=entry,                 
                description=input_dict[entry][0]                                                           
            )
            # Write the record to the output file in FASTA format
            SeqIO.write(record, output, "fasta")

def preparing_candidates_for_conflict_comparison(candidates_gff3, transcripts_gtf):
    chromosome_dict = {}
    with open(transcripts_gtf, "r") as transcripts_file:
        for line in transcripts_file:
            if line.startswith("#"):
                continue
            part = line.split('\t')
            if part[2] == "transcript":
                gene_id = re.search(r'gene_id "([^"]+)"', line)
                gene_id = gene_id.group(1)
                chromosome_dict[gene_id] = part[0]

    with open("candidates.bed", "w") as output, open(candidates_gff3, "r") as candidates_file:
        for line in candidates_file:
            #Bei Leerzeile, also wenn nicht strip() gemacht werden kann wird übersprungen
            if line.startswith("#") or not line.strip():
                continue
            part = line.split('\t')
            if part[2] == "CDS":
                start = part[3]
                stop = part[4]
                gene_id = re.search(r"(STRG\.\d+)", line).group(1)
                transcript_id = re.search(r'Parent=([^\s;]+)', line).group(1)
                chromosome = chromosome_dict[gene_id]
                strand = part[6]
                output.write(f"{chromosome}\t{start}\t{stop}\t{transcript_id}\t.\t{strand}\n")

def finding_protein_conflicts(candidates_bed, reference_bed):
    try:
        command = [
            "/home/s-amknut/GALBA/tools/bedtools2/bin/bedtools",
            "coverage",
            "-a",
            candidates_bed,
            "-b",
            reference_bed,
            "-s" #strand wird berücksichtigt
        ]
        with open("conflicts.bed", "w") as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE)

        if result.returncode == 0:
            print("Coverage file created successfully")

        else:
            print("Error during creating coverage file")
            print(result.stderr)

    except Exception:
        print("Could not run bedtools command.")
        sys.exit(1)

def finding_stop_in_utr(transcripts_fasta, intrinsic_candidates_genome):
    stop_in_utr_dict = {}   
    transcripts_dict = {}                   
    stop_codons = ["TAA", "TAG", "TGA"]
    for record in SeqIO.parse(transcripts_fasta, "fasta"):
        transcripts_dict[record.id] = [record.description, record.seq]
    with open(intrinsic_candidates_genome, "r") as candidates_file:
        for line in candidates_file:
            has_five_prime_utr = False
            if line.startswith("#"):
                continue

            part = line.strip().split("\t")
            if len(part) < 9:  # Skip empty lines
                continue
            feature_type = part[2]
            start = int(part[3])
            end = int(part[4])
            attributes = part[8]
            stringtie_id = re.search(r"ID=([^;]+?)\.p\d+\b", attributes)
            if stringtie_id:
                stringtie_id = stringtie_id.group(1)

            if feature_type == "five_prime_UTR":
                utr_length = end - start + 1
                utr_sequence = transcripts_dict[stringtie_id][1][utr_length:]
                for codon in stop_codons:                       
                    stop_codon_position = utr_sequence.find(codon)
                    if stop_codon_position != -1:
                        stop_in_utr_dict[stringtie_id] = True
                        break

    return stop_in_utr_dict 

def getting_hc_supported_by_proteins(diamond_tsv, transdecoder_pep, protein_file):
    t_length_dict = {}
    for record in SeqIO.parse(protein_file, "fasta"):
        t_length_dict[record.id] = len(record.seq) 
    q_dict = {}
    already_hc_genes = []
    count = 0
    for record in SeqIO.parse(transdecoder_pep, "fasta"):
        # t_length_dict[record.id] = len(record.seq) - 1 # -1 damit * nicht gezählt wird
        q_dict[record.id] = [record.description, record.seq] #46187 Elemente in t_dict
    print("Länge vom dict vor Filterung", len(q_dict))
    with open(diamond_tsv, "r") as tsv, open("hc_genes.pep", "w") as output:
        for line in tsv:    #7580 Elemente in tsv
            part = line.strip().split('\t')
            id = part[0] #Query ID
            protein_id = part[1] #Target ID
            aaident = float(part[2])
            align_length = int(part[3])
            #t_start = int(part[6]) #Bis 18.11. dachte ich es wäre t zuerst und dann q, dann auf website von diamond was anderes gelesen
            #t_end = int(part[7])
            #q_start = int(part[8])
            #q_end = int(part[9])
            q_start = int(part[6])
            q_end = int(part[7])
            t_start = int(part[8])
            t_end = int(part[9])
            #bitscore = int(part[11])
            if id in q_dict:
                #print("----ID: ", id, "Protein ID: ", protein_id, "Query Start: ", q_start, "Query End: ", q_end, "Target Start: ", t_start, "Target End: ", t_end ,"-----")
                record.id = id
                record.description = q_dict[id][0]
                record.seq = q_dict[id][1]
                #t_length = int(re.search(r"len:(\d+)", record.description).group(1)) #geht auch, weil len genau sequenzlänge der AS angibt
                q_length = len(record.seq) - 1 # -1 damit * nicht gezählt wird
                t_length = t_length_dict[protein_id]
                #print("Query Length: ", q_length, "Target Length: ", t_length)
                start_condition = (q_start - t_start)
                stop_condition = (q_length - q_end) - (t_length - t_end) 
                #print("Start Condition: ", start_condition, "Stop Condition: ", stop_condition)
                
                if "type:complete" in record.description: #28490 complete & hc von insgesamt 42797 complete candidates 
                    #Combi 0 und 0:    62.7     |    88.7    |
                    #Combi 3 und 10:   66.4     |    87.9    |
                    #Combi 6 und 21:   66.5     |    87.9    | 
                    #Combi 6 und 21 und t_start == 1:    64.3     |    88.2   
                    #Combi 6 und 21 und t_start == 1 und q_start == 1 und (t_length-t_end)==0:    61.8     |    88.8 
                    #Combi 6 und 21 und (t_length-t_end) == 0:     64.6     |    88.5  
                    #Combi 6 und 21 und t_start == 1 und q_start == 1 und (t_length-t_end)==0 und (q_length-q_end)==0:    61.6     |    88.8
                    #Combi 0 und 0 und (t_length-t_end)==0 und aaident > 99.5:    49.7     |    89.1    |
                    #Combi 6 und 21 und (t_length-t_end)==0 und aaident > 99.5:    52.7     |    89.1    |
                    #Combi 6 und 21 und aaident > 99.5:    54.4     |    88.9    |
                    #Combi 6 und 21 und aaident > 99.8:    41.6     |    88.8    |
                    #Combi 6 und 21 und aaident > 98:    62.7     |    88.5    |
                    #Combi 6 und 21 und aaident > 97:    63.7     |    88.5    |
                    #Combi 6 und 21 und aaident > 96.5:    64.0     |    88.4    |
                    #Combi 6 und 21 und aaident > 96:    64.3     |    88.4    
                    #Combi 6 und 21 und aaident > 95:    64.7     |    88.4    |
                    #Combi 6 und 21 und aaident > 95 und (t_length-t_end)<6:    62.9     |    88.8    |
                    #Combi 6 und 21 und aaident > 95 und (t_length-t_end)<11:   63.0     |    88.8    |
                    #Combi 6 und 21 und aaident > 95 und (t_length-t_end)<15:   63.1     |    88.8    |
                    #Combi 6 und 21 und aaident > 95 und (t_length-t_end)<18:   63.1     |    88.8    |
                    #Combi 6 und 21 und aaident > 95 und (t_length-t_end)<15 und t_start < 5:   60.3     |    88.9    |
                    #Combi 6 und 21 und aaident > 95 und (t_length-t_end)<15 und t_start < 25:     60.8     |    88.8    |
                    #Combi 0 und 0 und aaident == 100 und (t_length-t_end)==0 und q_start == 1:     31.6     |    88.1    
                    #and aaident > 95 and (t_length-t_end)<15
                    #if start_condition < 6 and stop_condition < 21: #abs() AUSTESTEN # Eigentlich 6 und 21 aber neu getestet mit 3 und 10
                    #if aaident > 95 and (t_length-t_end)<10:
                    if (t_length - align_length) < 15 and start_condition < 6 and aaident > 95:
                    #if start_condition < 6 and (t_length - align_length) < 15:
                        SeqIO.write(record, output, "fasta")
                        gene_id = id.split(".")[0] + "." + id.split(".")[1] 
                        del q_dict[id]
                        already_hc_genes.append(gene_id)
                if "type:5prime_partial" in record.description or "type:internal" in record.description or "type:3prime_partial" in record.description:
                    #print ("Nicht erfüllt wegen incomplete")
                    del q_dict[id]
                        #print("StringTie ID ist hc: ", stringtie_id)
                        #keys_to_delete = [key for key in t_dict if key.startswith(gene_id)]
                        #for key in keys_to_delete:
                            #print("Key wird gelöscht: ", key)
                           # del t_dict[key]
                        #continue
                    #else:
                     #   if gene_id not in hc_genes:   
                        #Hier intrinsic, complete schon abgehakt.
                #for id in q_dict:
                 #   gene_id = id.split(".")[0] + "." + id.split(".")[1] 
                  #  if gene_id in already_hc_genes:
                   #     print(gene_id)
                    #    del q_dict[id]
        print(len(q_dict))
        q_dict = {
            id: value
            for id, value in q_dict.items()
            if id.split(".")[0] + "." + id.split(".")[1] not in already_hc_genes
        }
        print(len(q_dict))

    return q_dict

def getting_hc_supported_by_intrinsic(q_dict):
    from_dict_to_pep_file(q_dict, "intrinsic_candidates.pep")
    choose_one_isoform("intrinsic_candidates.pep", "intrinsic_one_isoform.pep")
    from_pep_file_to_gff3("intrinsic_one_isoform.pep", "transcripts.gtf", "intrinsic_candidates.gff3")
    from_transcript_to_genome("intrinsic_candidates.gff3", "transcripts.gff3", "transcripts.fasta", "intrinsic_candidates_genome.gff3")
    preparing_candidates_for_conflict_comparison("intrinsic_candidates_genome.gff3", "transcripts.gtf")
    preparing_miniprot_gff_for_conflict_comparison("miniprot_parsed.gff")
    finding_protein_conflicts("candidates.bed", "reference.bed")
    stop_in_utr_dict= finding_stop_in_utr("transcripts.fasta", "intrinsic_candidates_genome.gff3")
    count2 = 0
    count1 = 0
    no_conflicts_dict = {}
    with open("conflicts.bed", "r") as conflicts_file:
        for line in conflicts_file:
            part = line.split('\t')
            transcript_id = part[3]
            if int(part[6]) == 0: 
                no_conflicts_dict[transcript_id] = True
            else:
                no_conflicts_dict[transcript_id] = False

    with open("hc_genes.pep", "a") as output, open("intrinsic_one_isoform.pep") as candidates_file:
        cond1_true = 0
        cond2_true = 0 
        cond3_true = 0
        cond4_true = 0 
        generell = 0
        bestehen = 0
        for record in SeqIO.parse(candidates_file, "fasta"):
            count1 += 1
            transcript_id = record.id
            length_pep = int(re.search(r"len:(\d+)", record.description).group(1))
            length_cds = length_pep*3
            score = float(re.search(r"score=(-?[\d.]+)", record.description).group(1))
            if transcript_id in no_conflicts_dict: #falls gegen Genom gemapped wurde
                condition1 = (length_cds >= 300)
                condition2 = (score > 50) 
                condition3 = (no_conflicts_dict[transcript_id])
                id_stringtie = transcript_id.split(".p")
                id_stringtie = id_stringtie[0]
                if id_stringtie in stop_in_utr_dict:
                    condition4 = (stop_in_utr_dict[id_stringtie])
                else:
                    condition4 = False
                if condition1 and condition2 and condition3 and condition4:
                    SeqIO.write(record, output, "fasta")
                    bestehen += 1
                if condition1 == True: 
                    cond1_true += 1
                if condition2 == True: 
                    cond2_true += 1
                if condition3 == True: 
                    cond3_true += 1
                if condition4 == True: 
                    cond4_true += 1
                generell += 1
    print("Score 70")
    print("Insgesamt: ", generell)
    print("Condition 1: ", cond1_true)
    print("Condition 2: ", cond2_true)
    print("Condition 3: ", cond3_true)
    print("Condition 4: ", cond4_true)
    print("HC: ", bestehen)

        
               
'''
        last_gene_id = None
        longest_length = 0
        intrisic_dict = {}
        for id in t_dict:
            record.id = id
            record.description = t_dict[id][0]
            record.seq = t_dict[id][1]
            length = int(re.search(r"len:(\d+)", record.description).group(1))
            gene_id = re.search(r"(STRG\.\d+)\.\d+", record.description).group(1)
            print("ID: ", id, " hat Länge: ", length)
            if gene_id == last_gene_id or last_gene_id == None:
                if length > longest_length:
                    print("Neues längstes Transkript: ", length, "ID: ", id)
                    longest_length = length
                    longest_id = id
            else:
                print("Längstes Transkript: ", longest_length, "ID: ", longest_id)
                intrisic_dict[longest_id] = t_dict[longest_id]
                longest_length = length
                longest_id = id
            last_gene_id = gene_id
'''
        
            
                #Für nicht complete cds: DRINGEND WIEDER HINZUFÜGEN 
               # if "type:5prime_partial" in record.description: #1021 5prime & hc von insgesamt 5126 5prime candidates
                #    if stop_condition < 21:
                 #       SeqIO.write(record, output, "fasta")
                  #      del t_dict[id]
                   #     continue

        #for id in t_dict:

def choose_one_isoform(transdecoder_pep, output_name): 
    isoform_dict = {}
    sum_of_scores = 0
    amount_of_entrys = 0
    smallest_score = 1000000
    with open(output_name, "w") as output:
        for record in SeqIO.parse(transdecoder_pep, "fasta"):
            transdecoder_id = record.id
            stringtie_id = transdecoder_id.split(".p")[0]
            gene_id = stringtie_id.split(".")[1]
            description = record.description.split(" ")
            cds_coords = re.search(r":(\d+)-(\d+)\((\+|\-)\)", description[7])
            start_cds = int(cds_coords.group(1))
            stop_cds = int(cds_coords.group(2))
            description = record.description.split(" ")

            score = re.search(r"score=(-?[\d.]+)", record.description)
            score = float(score.group(1))
            if gene_id not in isoform_dict:
                isoform_dict[gene_id] = [(start_cds, stop_cds, record, score)]
            else:
                isoform_dict[gene_id].append((start_cds, stop_cds, record, score))

        for gene_id in isoform_dict:
            if len(isoform_dict[gene_id]) == 1:
                record = isoform_dict[gene_id][0][2]
                SeqIO.write(record, output, "fasta")
                sum_of_scores += isoform_dict[gene_id][0][3]
                amount_of_entrys += 1
                if isoform_dict[gene_id][0][3] < smallest_score:
                    smallest_score = isoform_dict[gene_id][0][3]
            else:
                #takes the longest isoform, if there are two longest, it takes the first, because that indicates the highest expression rate
                longest_isoform = max(isoform_dict[gene_id], key=lambda x: x[1] - x[0]) 
                record = longest_isoform[2]
                SeqIO.write(record, output, "fasta")
                sum_of_scores += isoform_dict[gene_id][0][3]
                amount_of_entrys += 1
                ##orf_highest_score = max(isoform_dict[gene_id], key=lambda x: x[-1])
                ##record = orf_highest_score[2]
                ##SeqIO.write(record, output, "fasta")
    print("Sum scores: ", sum_of_scores) #1408841.4800000126
    print("Amount entrys: ", amount_of_entrys) #18844
    #Mean: 74.76
    print("Smallest Score: ", smallest_score)
    
def parse_transdecoder_file(transdecoder_pep):
    stringtie_id_dict = {} 
    for record in SeqIO.parse(transdecoder_pep, "fasta"):
        transdecoder_id = record.id
        stringtie_id = transdecoder_id.split(".p")[0]
        description = record.description.split(" ")
        cds_transcript_coords = re.search(r":(\d+)-(\d+)\((\+|\-)\)", description[7])
        start_cds_transcript = int(cds_transcript_coords.group(1))
        stop_cds_transcript = int(cds_transcript_coords.group(2))
        cds_length = int(stop_cds_transcript) - int(start_cds_transcript) + 1
        strand = re.search(r"\((\+|-)\)", description[7])
        strand = strand.group(1)
        triple = (start_cds_transcript, stop_cds_transcript, cds_length, strand) #TRIPLE UMBENNENEN
        #Hier noch einbauen, dass wir nur das längste nehmen, falls wir am Ende eine file haben wollen, die mehrere isoformen enthalten soll
        if stringtie_id not in stringtie_id_dict:
            stringtie_id_dict[stringtie_id] = [triple]
        else:
            stringtie_id_dict[stringtie_id].append(triple)
    filtered_dict = {stringtie_id: val for stringtie_id, val in stringtie_id_dict.items() if len(val) == 1} #nicht notwendig
    return filtered_dict

def parse_transdecoder_dict(transdecoder_dict):
    stringtie_id_dict = {}
    for record in transdecoder_dict:
        transdecoder_id = record
        #stringtie_id = transdecoder_id.split(".p")[0]
        description = transdecoder_dict[record][0]
        cds_transcript_coords = re.search(r":(\d+)-(\d+)\((\+|\-)\)", description)
        start_cds_transcript = int(cds_transcript_coords.group(1))
        stop_cds_transcript = int(cds_transcript_coords.group(2))
        cds_length = int(stop_cds_transcript) - int(start_cds_transcript) + 1
        triple = (start_cds_transcript, stop_cds_transcript, cds_length)
        if transdecoder_id not in stringtie_id_dict:
            stringtie_id_dict[transdecoder_id] = [triple]
        else:
            stringtie_id_dict[transdecoder_id].append(triple)
    return stringtie_id_dict

def from_transcript_to_genome(orf_gff3, transcripts_gff3, transcripts_fasta, output_name):
    try:
        command = [
            "/opt/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl",
            orf_gff3,
            transcripts_gff3,
            transcripts_fasta,
            #">",
            #output_name
        ]
        with open(output_name, "w") as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE)
       # result = subprocess.run(command)

        if result.returncode == 0:
            print("Successfull")

        else:
            print("Error")
            print(result.stderr)
    
    except Exception:
        print("Could not run cdna_alignment_orf_to_genome_orf.pl")
        sys.exit(1)

def from_transcript_to_genome_coords(stringtie_gtf, transdecoder_id_dict, output_file): #better name: creating_annotation_file
    count = 0
    with open(stringtie_gtf, "r") as stringtie, open(output_file, "w") as output:
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
                attributes = part[8]

                gene_id = re.search(r'gene_id "([^"]+)"', attributes)
                gene_id = gene_id.group(1)
                transcript_id = re.search(r'transcript_id "([^"]+)"', attributes)
                transcript_id = transcript_id.group(1)
                #Eventuell noch exon_number hinzufügen
                #Soll score mit rein? Habs erstmal rausgenommen

                if feature == 'transcript':
                    print("--------------Neues Transcript mit ID: ----------- ", transcript_id) 
                    output.write(f"{seqname}\tPreGalba\tgene\t{start_genome}\t{stop_genome}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    output.write(f"{seqname}\tPreGalba\ttranscript\t{start_genome}\t{stop_genome}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    #print("ANFANG: Transkriptlänge, cds Länge, cds index auf 0")
                    if len(exon_coords_list) > 1: #Für CDS muss es auch zugelassen sein, dass nur ein Exon vorhanden ist
                        transcript_id = exon_coords_list[0][2]
                        gene_id = exon_coords_list[0][3]
                        strand = exon_coords_list[0][4]
                        print("Transkript ID: ", transcript_id, "Exon Liste Länge: ", len(exon_coords_list), spliced_transcript_length)
                        #print("Exons für vorheriges Transkript vorhanden: ", transcript_id)
                        for i in range(len(exon_coords_list)-1):
                            curr_exon_stop = int(exon_coords_list[i][1])
                            next_exon_start = int(exon_coords_list[i+1][0])
                            intron_start = int(curr_exon_stop) + 1
                            intron_stop = int(next_exon_start) - 1
                            #print("Intron von: ", intron_start, "bis: ", intron_stop)
                            output.write(f"{seqname}\tPreGalba\tintron\t{intron_start}\t{intron_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    if len(exon_coords_list) > 0:
                        if transcript_id in transdecoder_id_dict:
                            strand_transcript = transdecoder_id_dict[transcript_id][0][3]
                            if strand_transcript == "+":
                                cds_transcript_start = int(transdecoder_id_dict[transcript_id][0][0])
                                cds_transcript_stop = int(transdecoder_id_dict[transcript_id][0][1])
                            if strand_transcript == "-":
                                print("Spliced Transcript: ", spliced_transcript_length)
                                print("Start in Transdecoder: ", transdecoder_id_dict[transcript_id][0][0])
                                print("Stop in Transdecoder: ", transdecoder_id_dict[transcript_id][0][1])
                                cds_transcript_stop = spliced_transcript_length - int(transdecoder_id_dict[transcript_id][0][0]) + 1
                                cds_transcript_start = spliced_transcript_length - int(transdecoder_id_dict[transcript_id][0][1]) + 1
                                print("Neuer Start: ", cds_transcript_start)
                                print("Neuer Stopp: ", cds_transcript_stop)
                            cds_total_length = int(transdecoder_id_dict[transcript_id][0][2])
                            
                           # print("CDS ist insgesamt so lang: ", cds_total_length)
                           # print("CDS Koordinaten in Transkript: ", cds_transcript_start, cds_transcript_stop)
                            curr_transcript_length = 0
                            cds_current_length = 0
                            for i in range(len(exon_coords_list)):
                                curr_exon_start = int(exon_coords_list[i][0]) 
                                curr_exon_stop = int(exon_coords_list[i][1])
                                #print("Exonkoordinaten vom aktuellen Exon: ", curr_exon_start, curr_exon_stop)
                                if i > 0:
                                    curr_transcript_length += int(exon_coords_list[i-1][1]) - int(exon_coords_list[i-1][0]) + 1
                                    #print("Vorherige Transkriptlänge wird hochgesetzt auf: ", curr_transcript_length)
                                #Anfang vom CDS:
                                if cds_current_length == 0:
                                    if curr_exon_start - curr_transcript_length + cds_transcript_start - 1 > curr_exon_stop: #GEÄNDERT
                                        x = curr_exon_start + cds_transcript_start - 1 #nochmal genau bestimmen GEÄNDERT
                                        #print("CDS Startpunkt:", x , " liegt hinter Exonstoppunkt: ", curr_exon_stop, "Also zu nächstem Exon springen")
                                        fivePrimeUTR_start = curr_exon_start
                                        fivePrimeUTR_stop = curr_exon_stop
                                        if strand == "+":
                                            output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        else:
                                            output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        continue
                                    else: 
                                        cds_start_genome = curr_exon_start + int(cds_transcript_start) - curr_transcript_length - 1
                                        #print("CDS Startpunkt liegt innerhalb des Exons bei: ", cds_start_genome)
                                        if cds_start_genome + cds_total_length - cds_current_length - 1 > curr_exon_stop: #GEÄNDERT
                                            #print("CDS geht bis Exongrenze...")
                                            output.write(f"{seqname}\tPreGalba\tCDS\t{cds_start_genome}\t{curr_exon_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            #frame = frame(frame, curr_exon_stop - cds_start_genome + 1, strand)
                                            #print("CDS von: ", cds_start_genome, "bis: ", curr_exon_stop, "hinzugefügt")
                                            cds_current_length = curr_exon_stop - cds_start_genome + 1
                                            #print("CDS aktuelle Länge: ", cds_current_length)
                                            fivePrimeUTR_start = curr_exon_start
                                            fivePrimeUTR_stop = cds_start_genome - 1
                                            start_codon_plus_start = cds_start_genome
                                            start_codon_plus_stop = cds_start_genome + 2
                                            if strand == "+":
                                                output.write(f"{seqname}\tPreGalba\tstart_codon\t{start_codon_plus_start}\t{start_codon_plus_stop}\t.\t{strand}\t0\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                if curr_exon_start != cds_start_genome:
                                                    output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            else:
                                                output.write(f"{seqname}\tPreGalba\tstop_codon\t{start_codon_plus_start}\t{start_codon_plus_stop}\t.\t{strand}\t0\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                if curr_exon_start != cds_start_genome:
                                                    output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")

                                        else:
                                            #print("CDS innerhalb des Exons beendet...")
                                            cds_stop_genome = cds_start_genome + cds_total_length - 1
                                            cds_current_length = cds_total_length
                                            #print("Aktueller Exonstoppunkt: ", stop_genome)
                                            #print("CDS Stoppunkt: ", cds_stop_genome)
                                            output.write(f"{seqname}\tPreGalba\tCDS\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            fivePrimeUTR_start = curr_exon_start
                                            fivePrimeUTR_stop = cds_start_genome - 1
                                            threePrimeUTR_start = cds_stop_genome + 1
                                            threePrimeUTR_stop = curr_exon_stop
                                            #print("CDS von: ", cds_start_genome, "bis: ", cds_stop_genome, "hinzugefügt")
                                            start_codon_plus_start = cds_start_genome
                                            start_codon_plus_stop = cds_start_genome + 2
                                            stop_codon_plus_start = cds_stop_genome - 2
                                            stop_codon_plus_stop = cds_stop_genome
                                            if strand == "+":
                                                if curr_exon_start != cds_start_genome:
                                                    output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                if curr_exon_stop != cds_stop_genome:
                                                    output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{threePrimeUTR_start}\t{threePrimeUTR_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                output.write(f"{seqname}\tPreGalba\tstart_codon\t{start_codon_plus_start}\t{start_codon_plus_stop}\t.\t{strand}\t0\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                output.write(f"{seqname}\tPreGalba\tstop_codon\t{stop_codon_plus_start}\t{stop_codon_plus_stop}\t.\t{strand}\t0\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            else:
                                                if curr_exon_start != cds_start_genome:
                                                    output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{fivePrimeUTR_start}\t{fivePrimeUTR_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                if curr_exon_stop != cds_stop_genome:
                                                    output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{threePrimeUTR_start}\t{threePrimeUTR_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                output.write(f"{seqname}\tPreGalba\tstop_codon\t{start_codon_plus_start}\t{start_codon_plus_stop}\t.\t{strand}\t0\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                                output.write(f"{seqname}\tPreGalba\tstart_codon\t{stop_codon_plus_start}\t{stop_codon_plus_stop}\t.\t{strand}\t0\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            continue
                                        
                                elif cds_current_length == cds_total_length:
                                    #print("CDS bereits beendet. Rest der exons wird zu UTR")
                                    if strand == "+":
                                        output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{curr_exon_start}\t{curr_exon_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                    else:
                                        output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{curr_exon_start}\t{curr_exon_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                else:
                                    cds_start_genome = curr_exon_start
                                    if cds_current_length+int(curr_exon_stop)-int(curr_exon_start)+1<cds_total_length:
                                        #print("CDS noch nicht beendet...")
                                        cds_stop_genome = int(curr_exon_stop) 
                                        #print("CDS wird von Start des aktuellen Exons, bis Stopp des aktuellen Exons eingetragen")
                                        output.write(f"{seqname}\tPreGalba\tCDS\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        cds_current_length+=curr_exon_stop - curr_exon_start + 1
                                        #print("Neue aktuelle Länge: ", cds_current_length)
                                    else:
                                        #print("CDS innerhalb des Exons beendet...")
                                        cds_stop_genome = int(curr_exon_start) + cds_total_length - cds_current_length - 1 #GEÄNDERT (-1 hinzugefügt)
                                        #print("Aktueller Exonstoppunkt: ", curr_exon_stop)
                                        #print("CDS Stoppunkt: ", cds_stop_genome)
                                        output.write(f"{seqname}\tPreGalba\tCDS\t{cds_start_genome}\t{cds_stop_genome}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        threePrimeUTR_start = cds_stop_genome + 1
                                        threePrimeUTR_stop = curr_exon_stop
                                        stop_codon_plus_start = cds_stop_genome - 2
                                        stop_codon_plus_stop = cds_stop_genome
                                        if strand == "+":
                                            output.write(f"{seqname}\tPreGalba\tstop_codon\t{stop_codon_plus_start}\t{stop_codon_plus_stop}\t.\t{strand}\t0\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            if curr_exon_stop != cds_stop_genome:
                                                output.write(f"{seqname}\tPreGalba\tthree_prime_UTR\t{threePrimeUTR_start}\t{threePrimeUTR_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        else:
                                            output.write(f"{seqname}\tPreGalba\tstart_codon\t{stop_codon_plus_start}\t{stop_codon_plus_stop}\t.\t{strand}\t0\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                            if curr_exon_stop != cds_stop_genome:
                                                output.write(f"{seqname}\tPreGalba\tfive_prime_UTR\t{threePrimeUTR_start}\t{threePrimeUTR_stop}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                                        continue
                                        
                    exon_coords_list.clear()
                    spliced_transcript_length = 0
                        
                if feature == 'exon':
                    output.write(f"{seqname}\tPreGalba\texon\t{start_genome}\t{stop_genome}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";\n")
                    exon_number = re.search(r'exon_number "([^"]+)"', attributes)
                    exon_number = exon_number.group(1)
                    exon_coords_list.append((start_genome, stop_genome, transcript_id, gene_id, strand))
                    print(transcript_id, "Exon: ", start_genome, stop_genome)
                    spliced_transcript_length += (int(stop_genome) - int(start_genome)) + 1
                    count += 1

def frame_in_annotation(annotation_file):
    try:
        command = [
            "/home/s-amknut/GALBA/tools/eviann/src/gffread", #Schauen ob in sif file
            "-F",
            annotation_file,
            "-o",
            "annotation_with_frame.gff3"
        ]
        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Frames were successfully added to annotation file.")
        else:
            print("Error during adding frames to annotation file with gffread.")
            print(result.stderr)
    
    except Exception:
        print("Could not run gffread command.")
        sys.exit(1)

def only_cds_in_annotation(annotation_file):
    try:
        command = "grep 'CDS' " + annotation_file + " > annotation_only_cds.gff3" 

        result = os.system(command)
    
    except Exception:
        print("Could not run grep command.")
        sys.exit(1)
    
def control_annotation(annotation_file, reference, projname):
    try:
        command = [
            "/home/s-amknut/GALBA/tools/gffcompare-0.12.6.Linux_x86_64/gffcompare",
            "-r",
            reference,
            "-o",
            projname,
            annotation_file
        ]
    
        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Annotation stats were determined successfully.")
        else:
            print("Error during determining annotation stats with gffcompare.")
            print(result.stderr)

    except Exception:
        print("Could not run gffcompare command.")
        sys.exit(1)

def load_config(config_file):
    with open(config_file, "r") as config_file:
        input_files = yaml.safe_load(config_file)
        return input_files

#MAIN
parser = argparse.ArgumentParser(description='Genome annotation with transcriptomic data like RNA-seq and Iso-seq data')  
parser.add_argument('-t', '--threads', default=4, help='Number of threads (default=4)', required=False)
parser.add_argument('-y', '--config', help='Config file input', metavar='<config.yaml>', required=False) #required=True
parser.add_argument('-g', '--genome', help='Genome file input', metavar='<genome.fasta>', required=False)
parser.add_argument('-p', '--proteins', help='Protein file input', metavar='<proteins.fasta>', required=False)

parser.add_argument('--isoseq', action='store_true', help='Use this option if you want to process isoseq data only')
parser.add_argument('--rnaseq', action='store_true', help='Use this option if you want to process rnaseq data only')
parser.add_argument('--mixed', action='store_true', help='Use this option if you want to process both rnaseq and isoseq data')
parser.add_argument('--projname', default="pregalba", help='Name the output folder', required=False)
parser.add_argument('--output_path',  help='Specify the path you want the output folder to be created in if it should differ from the current working directory', required=False)

args = parser.parse_args()
threads = args.threads
proj_name = args.projname
output_path = args.output_path
print("Projektname: ", proj_name)
print("Parameter für HC durchtesten.")
if output_path == None:
    output_path = os.getcwd()

os.makedirs(output_path + "/" + proj_name, exist_ok=True)
os.chdir(output_path + "/" + proj_name)

if args.config:
    input_files = load_config(args.config)
    genome_file = input_files["genome"]
    rnaseq_paired_sets = input_files.get("rnaseq_paired_sets", []) #Wenn Liste nicht vorhanden, dann leere Liste
    rnaseq_single_sets = input_files.get("rnaseq_single_sets", [])
    isoseq_sets = input_files.get("isoseq_sets", [])
    protein_file = input_files["protein"] #optional?
    reference_annotation = input_files["annotation"]

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

if args.genome:
    genome_file = args.genome

if args.proteins:
    protein_file = args.proteins
'''
if process_rnaseq:
    indexing(genome_file)
    alignments_list = mapping_short(rnaseq_paired_sets, rnaseq_single_sets)
    sam_to_bam(alignments_list)    
    if len(alignments_list) > 1:    #NOCHMAL PRÜFEN OB HIER WIRKLICH BAMFILES ÜBERGEBEN WERDEN
        alignment_rnaseq = file_name(merge_bam_files(alignments_list[0], alignments_list[1])) + ".bam"
    else:
        alignment_rnaseq = file_name(alignments_list[0]) + ".bam"
else:
    alignment_rnaseq = None

if process_isoseq:
    alignment_isoseq = mapping_long(genome_file, isoseq_sets)
    sam_file_list = [alignment_isoseq]
    sam_to_bam(sam_file_list) 
    alignment_isoseq = file_name(alignment_isoseq) + ".bam"
else:
    alignment_isoseq = None
'''
#assembling(alignment_rnaseq, alignment_isoseq)  #Für alleine testen leer machen
#orfsearching(genome_file, "transcripts.gtf")  #Vielleicht eher Was returned wurde als input übergeben
#convert_gtf_to_gff3("transcripts.gtf")
#short_start_dict = shorten_incomplete_Orfs("transcripts.fasta.transdecoder.pep")
#make_diamond_db(protein_file)
#validating_ORFs("shortened_candidates.pep", "diamond_shortened.tsv")
#make_diamond_db(protein_file)
#validating_ORFs("transcripts.fasta.transdecoder.pep", "diamond_normal.tsv")
#classifications_dict = get_cds_classification("diamond_normal.tsv", "diamond_shortened.tsv", short_start_dict)
#get_optimized_pep_file("transcripts.fasta.transdecoder.pep", "shortened_candidates.pep", classifications_dict, short_start_dict)
#make_diamond_db(protein_file)
#validating_ORFs("revised_candidates.pep", "diamond_revised.tsv")

q_dict = getting_hc_supported_by_proteins("diamond_revised.tsv", "revised_candidates.pep", protein_file)
#protein_aligning(genome_file, protein_file, "/home/s-amknut/GALBA/tools/blosum62_1.csv") 
getting_hc_supported_by_intrinsic(q_dict)

choose_one_isoform("hc_genes.pep", "one_chosen_isoform.pep")
#choose_one_isoform("revised_candidates.pep", "one_chosen_isoform.pep")
#choose_one_isoform("transcripts.fasta.transdecoder.pep", "one_chosen_isoform.pep")
#transdecoder_id_dict = parse_transdecoder_file("one_chosen_isoform.pep")
#from_transcript_to_genome_coords("transcripts.gtf", transdecoder_id_dict, "annotation.gtf")
from_pep_file_to_gff3("one_chosen_isoform.pep", "transcripts.gtf", "one_chosen_isoform.gff3")
#from_transcript_to_genome("transcripts.fasta.transdecoder.gff3","transcripts.gff3","transcripts.fasta", "transcripts.fasta.transdecoder.genome.gff3")
from_transcript_to_genome("one_chosen_isoform.gff3","transcripts.gff3","transcripts.fasta", "transcripts.fasta.transdecoder.genome.gff3")
###frame_in_annotation("transcripts.fasta.transdecoder.genome.gff3")
#only_cds_in_annotation("annotation_with_frame.gff3")
only_cds_in_annotation("transcripts.fasta.transdecoder.genome.gff3")
control_annotation("annotation_only_cds.gff3", reference_annotation, proj_name) #noch testen
#frame_in_annotation("annotation.gtf")
#only_cds_in_annotation("annotation_with_frame.gtf")
#control_annotation("annotation_only_cds.gtf", reference_annotation, proj_name) #noch testen
#only_cds_in_annotation("transcripts.fasta.transdecoder.genome.gff3")
#control_annotation("annotation_only_cds.gtf", reference_annotation, proj_name) #noch testen
#control_annotation("/home/s-amknut/GALBA/braker/GeneMark-ETP/training_cds.gtf", reference_annotation, "new_braker_cds")
#preparing_miniprot_gff_for_conflict_comparison("miniprot_parsed.gff")
#preparing_candidates_for_conflict_comparison(t_dict)
#finding_protein_conflicts("bedtools_reference.bed", "bedtools_candidates.bed")



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
#-Sind in der output file noch gene ohne cds? Wenn ja rausnehmen

#FRAGEN:
#-Sollte ich mit -G die stringtie Option nutzen, eine Referenzannotation zu verwenden? --> Diese dann in die yaml file oder parser?
#-Optionen in Ordnung oder noch Unterscheidung zwischen single und paired-end?
#-Ist es richtig, dass man single-end und paired-end beide nutzt? Oder wird in der Regel nur eins davon genutzt?
#-Soll die Übergabe von Proteinfiles optional sein?
#-Richtig, dass alle erstellten files in cwd gespeichert werden? Soll ich Funktion einfügen, dass man sich das aussuchen kann wohin?
#-Scoring Matrix von Nutzer einfügen lassen oder selbst eine vorgeben?
#-Ist es richtig, dass ein Trankript nach dem splicing mehrere CDS haben kann? D.h. dass man UTR mittig hat?

'''Neue Fragen'''
#-Wie bekomme ich von der miniprot file die Chromosomennummer
#-Bedtools okay?
#Betrag bei candidates 
#-Wähle ich zu viele HC Gene aus #28490 complete & hc von insgesamt 42797 complete candidates 
#-Gibt es eine Möglichkeit, dass ich einzelne mit dem von GeneMark vergleichen kann?

#Plan:
#-Transdecoder macht ORF prediction -> Davor intron hints von spliced rnaseq daten wie bei genemark?
#-miniprot macht protein alignment
#-Proteinalignment von miniprot mit dem von Transdecoder vergleichen -> Prediction ergänzen oder verwerfen???
#-Mit miniprot trainingsgenen Augustus trainieren und hc hints an Augustus übergeben

#Oder:
#-Transdecoder macht ORF prediction -> Davor intron hints von spliced rnaseq daten wie bei genemark?
#-Diamond nutzt .pep file von Transdecoder und sucht homologe Proteine 
#-Spaln aligniert die homologen Proteine zurück ans Genom, um predictions genauer zu machen 

#File, die training.gtf entspricht + file mit allen hc genen + file mit allen genen + file mit allen hc introns
#Nur eine Isoform pro Gen wählen: Längste CDS und sonst niedrigste StringTie Nummer 
#CDS ohne Protein support 

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

training.gtf mit allen rnaseq und protein.fa Daten:
#-----------------| Sensitivity | Precision  |
        Base level:    66.3     |    99.5    |
        Exon level:    64.6     |    99.1    |
      Intron level:    70.3     |    99.4    |
Intron chain level:    40.1     |    96.4    |
  Transcript level:    39.5     |    96.2    |
       Locus level:    58.0     |    96.3    |

#training_cds.gtf mit allen rnaseq und protein.fa Daten:
#-----------------| Sensitivity | Precision  |
        Base level:    66.3     |    99.5    |
        Exon level:    64.6     |    99.1    |
      Intron level:    70.3     |    99.4    |
Intron chain level:    40.1     |    96.4    |
  Transcript level:    39.5     |    96.2    |
       Locus level:    58.0     |    96.3    |

#main.py mit allen rnaseq und isoseq Daten (mit Kategorisierung und mit hc Filter(ohne intrinsic und incomplete)), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Korrektur vom 17.11.:     
#-----------------| Sensitivity | Precision  |
        Base level:    63.9     |    88.5    |
        Exon level:    60.0     |    84.5    |
      Intron level:    72.3     |    91.9    |
Intron chain level:    27.8     |    62.4    |
  Transcript level:    24.3     |    61.1    |
       Locus level:    35.7     |    61.2    |

#main.py mit allen rnaseq und isoseq Daten (mit Kategorisierung und mit hc Filter(ohne intrinsic und incomplete, ohne abs())), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Korrektur vom 17.11.:
#-----------------| Sensitivity | Precision  |
        Base level:    65.1     |    88.3    |
        Exon level:    60.9     |    84.3    |
      Intron level:    73.5     |    91.7    |
Intron chain level:    28.1     |    62.0    |
  Transcript level:    24.6     |    60.6    |
       Locus level:    36.1     |    60.7    |

#main.py mit allen rnaseq und isoseq Daten (OHNE Kategorisierung und mit hc Filter(ohne intrinsic und incomplete, mit abs())), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Korrektur vom 17.11.:
#-----------------| Sensitivity | Precision  |
        Base level:    59.9     |    88.7    |
        Exon level:    56.1     |    84.8    |
      Intron level:    67.6     |    92.0    |
Intron chain level:    25.6     |    62.2    |
  Transcript level:    22.4     |    60.9    |
       Locus level:    32.9     |    61.0    |

#main.py mit allen rnaseq und isoseq Daten (OHNE Kategorisierung und mit hc Filter(ohne intrinsic aber mit ALLEN incomplete, mit abs())), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Korrektur vom 17.11.:
#-----------------| Sensitivity | Precision  |
        Base level:    65.8     |    87.5    |
        Exon level:    61.2     |    83.2    |
      Intron level:    74.2     |    91.2    |
Intron chain level:    28.5     |    61.2    |
  Transcript level:    24.9     |    59.6    |
       Locus level:    36.6     |    59.7    |

#main.py mit allen rnaseq und isoseq Daten (OHNE Kategorisierung und OHNE hc Filter), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Korrektur vom 17.11.:
#-----------------| Sensitivity | Precision  |
        Base level:    70.6     |    87.0    |
        Exon level:    64.7     |    82.4    |
      Intron level:    78.7     |    91.1    |
Intron chain level:    30.7     |    60.3    |
  Transcript level:    26.8     |    58.4    |
       Locus level:    39.3     |    58.5    |

#main.py mit allen rnaseq und isoseq Daten (mit Kategorisierung und OHNE hc Filter(ohne intrinsic und incomplete)), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Korrektur vom 17.11.:
#-----------------| Sensitivity | Precision  |
        Base level:    70.6     |    87.4    |
        Exon level:    64.9     |    82.8    |
      Intron level:    78.7     |    91.2    |
Intron chain level:    30.8     |    60.5    |
  Transcript level:    26.9     |    58.6    |
       Locus level:    39.4     |    58.7    |
    
#main.py mit allen rnaseq und isoseq Daten (mit Kategorisierung und OHNE hc Filter(ohne intrinsic und incomplete)), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Diamond-Spalten Korrektur vom 18.11.:
#-----------------| Sensitivity | Precision  |
        Base level:    70.6     |    87.3    |
        Exon level:    64.9     |    82.7    |
      Intron level:    78.7     |    91.2    |
Intron chain level:    30.8     |    60.4    |
  Transcript level:    26.8     |    58.5    |
       Locus level:    39.4     |    58.6    |

#main.py mit allen rnaseq und isoseq Daten (mit Kategorisierung und OHNE hc Filter(ohne intrinsic und incomplete)), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit + short_start Korrektur und nur cds vom 27.11.:
#-----------------| Sensitivity | Precision  |
        Base level:    70.6     |    87.4    |
        Exon level:    65.0     |    82.8    |
      Intron level:    78.7     |    91.2    |
Intron chain level:    30.9     |    60.6    |
  Transcript level:    26.9     |    58.7    |
       Locus level:    39.5     |    58.9    |
       
#main.py mit allen rnaseq und isoseq Daten (mit Kategorisierung und mit hc Filter(ohne intrinsic und incomplete)), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Diamond-Spalten Korrektur vom 18.11.: 
#-----------------| Sensitivity | Precision  |
        Base level:    65.7     |    88.0    |
        Exon level:    61.3     |    83.7    |
      Intron level:    74.0     |    91.2    |
Intron chain level:    27.8     |    60.5    |
  Transcript level:    24.3     |    58.8    |
       Locus level:    35.7     |    59.0    |

#main.py mit allen rnaseq und isoseq Daten (mit Kategorisierung und mit hc Filter(ohne intrinsic und incomplete)), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Diamond-Spalten Korrektur + short Korrektur 27.11.: 
#-----------------| Sensitivity | Precision  |
        Base level:    66.5     |    87.9    |
        Exon level:    62.0     |    83.6    |
      Intron level:    74.8     |    91.2    |
Intron chain level:    28.3     |    60.7    |
  Transcript level:    24.8     |    59.0    |
       Locus level:    36.4     |    59.2    |

#main.py mit allen rnaseq und isoseq Daten (mit Kategorisierung und mit hc Filter(ohne intrinsic und incomplete)), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit (-) STrang in Annotation Korrektur 28.11.: 
#-----------------| Sensitivity | Precision  |
        Base level:    69.5     |    90.3    |
        Exon level:    69.4     |    93.5    |
      Intron level:    77.2     |    93.9    |
Intron chain level:    33.0     |    70.5    |
  Transcript level:    29.0     |    69.0    |
       Locus level:    42.6     |    69.5    |

#main.py mit allen rnaseq und isoseq Daten (mit Kategorisierung und mit hc Filter mit abs()(ohne intrinsic und incomplete)), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Diamond-Spalten Korrektur + short Korrektur + cds vom 27.11.:
#-----------------| Sensitivity | Precision  |
        Base level:    63.9     |    88.5    |
        Exon level:    59.9     |    84.5    |
      Intron level:    72.2     |    91.8    |
Intron chain level:    27.9     |    62.5    |
  Transcript level:    24.4     |    61.2    |
       Locus level:    35.8     |    61.3    |

#main.py mit rnaseq (mit Kategorisierung und mit hc Filter(ohne intrinsic und incomplete)), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Korrektur vom 17.11.:
#-----------------| Sensitivity | Precision  |
        Base level:    57.8     |    89.6    |
        Exon level:    55.3     |    85.2    |
      Intron level:    66.4     |    92.0    |
Intron chain level:    25.4     |    63.6    |
  Transcript level:    22.2     |    62.2    |
       Locus level:    32.7     |    62.3    |

#main.py mit rnaseq (mit Kategorisierung und mit hc Filter(ohne intrinsic und incomplete)), nur CDS in file, nur eine Isoform pro Gen mit FRAME mit Korrektur vom 28.11.:     
#-----------------| Sensitivity | Precision  |
        Base level:    63.0     |    91.4    |
        Exon level:    63.5     |    93.4    |
      Intron level:    71.0     |    94.0    |
Intron chain level:    29.9     |    71.1    |
  Transcript level:    26.3     |    69.5    |
       Locus level:    38.6     |    69.9    |

#rnaseq combi 0 und 0
#-----------------| Sensitivity | Precision  |
        Base level:    59.0     |    92.1    |
        Exon level:    60.3     |    95.1    |
      Intron level:    67.1     |    95.1    |
Intron chain level:    28.9     |    74.4    |
  Transcript level:    25.4     |    73.2    |
       Locus level:    37.3     |    73.5    |

#rnaseq Combi 6 und 21 und aaident > 95 und (t_length-t_end)<15:
#-----------------| Sensitivity | Precision  |
        Base level:    60.1     |    92.3    |
        Exon level:    61.6     |    94.4    |
      Intron level:    68.9     |    94.7    |
Intron chain level:    29.1     |    73.1    |
  Transcript level:    25.6     |    71.8    |
       Locus level:    37.5     |    72.1    |

#rnaseq Combi 6 und 21 und (t_length-t_end)<15:
#-----------------| Sensitivity | Precision  |
        Base level:    61.6     |    91.8    |
        Exon level:    62.7     |    94.0    |
      Intron level:    70.0     |    94.4    |
Intron chain level:    29.6     |    72.2    |
  Transcript level:    26.0     |    70.8    |
       Locus level:    38.2     |    71.1    |

#mixed mit transdecoder genome coords funktion und nur transdecoder outout ohne categories und hc alle isoforms
#-----------------| Sensitivity | Precision  |
        Base level:    81.3     |    96.0    |
        Exon level:    79.0     |    87.3    |
      Intron level:    85.6     |    94.3    |
Intron chain level:    57.7     |    50.0    |
  Transcript level:    55.2     |    49.7    |
       Locus level:    70.4     |    89.5    |

#mixed mit transdecoder genome coords funktion und nur transdecoder mit categories und ohne hc, nur eine isoform:
#-----------------| Sensitivity | Precision  |
        Base level:    78.3     |    97.7    |
        Exon level:    74.0     |    97.0    |
      Intron level:    81.1     |    99.0    |
Intron chain level:    46.1     |    91.8    |
  Transcript level:    45.1     |    89.5    |
       Locus level:    66.3     |    89.6    |

#mixed mit transdecoder genome coords funktion und nur transdecoder mit categories und hc, nur eine isoform:
#-----------------| Sensitivity | Precision  |
        Base level:    73.7     |    98.2    |
        Exon level:    70.6     |    97.9    |
      Intron level:    77.0     |    99.0    |
Intron chain level:    42.8     |    92.6    |
  Transcript level:    42.1     |    90.8    |
       Locus level:    61.9     |    90.9    |

#mixed mit transdecoder genome coords funktion und nur transdecoder mit categories und hc (protein & intrinsic mit cond1 und 3), nur eine isoform:
#-----------------| Sensitivity | Precision  |
        Base level:    73.8     |    98.2    |
        Exon level:    70.6     |    97.9    |
      Intron level:    77.1     |    99.0    |
Intron chain level:    42.9     |    92.6    |
  Transcript level:    42.2     |    90.7    |
       Locus level:    61.9     |    90.8    |
    
#mixed mit transdecoder genome coords funktion und nur transdecoder mit categories und hc (protein & intrinsic mit cond1, cond2=70 und cond3), nur eine isoform:
#-----------------| Sensitivity | Precision  |
        Base level:    73.8     |    98.2    |
        Exon level:    70.6     |    97.9    |
      Intron level:    77.1     |    99.0    |
Intron chain level:    42.9     |    92.6    |
  Transcript level:    42.2     |    90.7    |
       Locus level:    61.9     |    90.8    |

#mixed mit transdecoder genome coords funktion und nur transdecoder mit categories und hc (start< 6 and stop< 21 aaident > 95 (t_length-t_end)<15 t_start < 5:), nur eine isoform:     
#-----------------| Sensitivity | Precision  |
        Base level:    66.7     |    99.4    |
        Exon level:    65.5     |    99.2    |
      Intron level:    71.4     |    99.5    |
Intron chain level:    40.5     |    96.7    |
  Transcript level:    39.5     |    96.2    |
       Locus level:    58.1     |    96.2    |

-Gibt hier auch unter hc noch ein paar wenige doppelte cds pro gen. Die werden rausgelöscht aus output file
    -> Idee: Wenn zwei hc cds in einem Gen sind, dann die mit dem besseren score behalten
-in trainigsgene file normalerweise (transcript, gene, exon, cds) und in hints (introns, start, stop). Ich habe jetzt 
    gerade alles in einem. Brauche ich auch beides?
-Gar kein Unterschied erkennbar bei Kategorisierung oder keiner Kategorisierung
-Wie kann es sein, dass die Exon level so viel schlechter sind und sollten wir sensitivität erhöhen?
-HC Auswahl bezieht sich gerade ausschließlich auf die CDS, weil wir dafür Proteindaten nutzen.
    Schmeißen wir transkript/gen raus wenn keine cds gefunden wurden oder sie nicht als hc kategorisiert wurden?
-Score für nicht Proteingestützte CDS: 
    -Inframe Stopcodon in 5' UTR verstehe ich nicht weil 1. Stopcodon in CDS und 2. Stopcodon in 3' UTR
    -GMS-T log-odds score > 50, ersetzen durch bit score? Z.B. Median berechnen von den anderen HC Genen und das festlegen als Grenze
    -Plan:
        -Habe hc.gff file von miniprot, hier stehen introns/start/stop drin
        -Nehme mir die exon/intron Struktur von Transkript wo CDS drin ist, ohne CDS ansich zu berücksichtigen
        -Suche danach wo mein Kandidat liegen könnte. 
        -Wenn es mehrere passende abschnitte gibt, nehme ich Protein mit dem besten score.
        -Wenn es keine Konflikte gibt ist es approved 

NEU
-Wie kann es sein, dass in Stringtie file die verschiedenen Isoformen leicht unterschiedliche Transkriptgrenzen haben?
-Entscheiden uns für Isoform mit längster CDS. Dadurch gibt tomas ja quasi auch vor, dass wir uns bei mehreren CDS im selben Transkript immer für die längste entscheiden.
Schließlich entscheiden wir uns ja nur wegen dieser CDS für die Isoform.
-Frames hinterfragen.
-Wozu TSEBRA, wenn wir schon eigentlich die Isoformen schon von Stringtie bekommen?
-Log-odds score 
-Abgabedatum passt mit Korrektur?
-Finale files: 
    -gene_models: Alle HC Genstrukturen
    -hints_file: Wie die von Tomas, nur CDS,Introns,Start,Stop (Keine Incomplete), habe ja jetzt gerade nur CDS, ändern?
    -Introns: Alle HC Introns
    -Alle HC und Low Confidence CDS
    Wo sind die incomplete HC? 
    Gibt es file mit allen isoforms
-CDS getrimmt?
-Kategorisierung macht gerade keinen Unterschied 

27.11.:
-Was kann ich zu den Daten schreiben?
-Was zitiere ich bei BRAKER3
-Kommt Pandas und Biopython mit in Tabelle?
-Dockerfile in Material erklären?
-Graphic tool / VS code und so erwähnen?
-Wie noch confidence erhöhen?


Zu tun:
Bedtools in die Dockerfile schreiben
miniprot in dockerfile nutzen
'''