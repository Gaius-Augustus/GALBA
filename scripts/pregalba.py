#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import re
import yaml
import pandas as pd
import math
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#Author: "Amrei Knuth"
#Credits: "Katharina Hoff"
#Version: "1.0"
#Email:"amrei.knuth@gmail.com"
#Date: "Janurary 2025"

#USAGE with docker container: 
#First pull container: singularity build galba.sif docker://amreiknuth/galba-notebook:devel 
#Second run pregalba.py: singularity exec -B ${PWD}:${PWD} galba.sif pregalba.py -y config.yaml

''' Getting the file format of a FASTA or FASTQ input file. '''
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

''' Getting the file name without the file extension. '''
def file_name(path): 
    file = path.split("/")[-1]
    name = file.split(".")[0]
    if name.endswith("_1") or name.endswith("_2"):
        prefix = name.split("_")[0]
        return prefix
    else:
        return name

''' Indexing with hisat2-build in preparation of the short read mapping '''
def indexing(genome_fasta):
    try:
        path = hisat + "/hisat2-build"
        command = [
            path, 
            "--quiet",
            "-p",
            str(threads),        
            genome_fasta,           
            "genome"   
        ]

        print("Building genome index...")
    
        result = subprocess.run(command, capture_output=False)
    
        if result.returncode == 0:
            print("Indexing completed successfully.")
        else:
            print("Error during indexing:")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)

    except Exception:
        print("Could not run hisat2-build command.") 
        sys.exit(1)
        
'''Mapping short reads to the genome using Hisat2'''
def mapping_short(rnaseq_paired_sets, rnaseq_single_sets):
    alignments_list = []
    path = hisat + "/hisat2" 
    try:
        if not rnaseq_single_sets == []:
            #Selecting an option for the hisat2 command based on the file format of the input files
            if file_format(rnaseq_single_sets[0]) == "fasta": 
                format_option = "-f"
            if file_format(rnaseq_single_sets[0]) == "fastq":
                format_option = ""
            if file_format(rnaseq_single_sets[0]) == "unknown":
                print("Error: Unknown file format for single-end short reads. Please provide a FASTA or FASTQ file.")
                sys.exit(1)
            string_with_sets = ",".join(rnaseq_single_sets)
            output_1 = "alignment_single_rnaseq.sam" 
            hisat2_single_command = [
                    path,  
                    format_option,                       
                    "-x", "genome",            
                    "-U", string_with_sets,
                    "--dta",
                    "-p", str(threads),
                    "-S", output_1,                  
                ]
            print("Mapping single-end short reads to the genome...")
            
            result = subprocess.run(hisat2_single_command, capture_output=True)
            
            if result.returncode == 0:
                print("Mapping of single-end short reads completed successfully.")
                alignments_list.append(output_1)

            else:
                print("Error during mapping of single-end short-reads:")
                print(result.stdout)
                print(result.stderr)
                sys.exit(1)
        
    except Exception:
        print("Could not run hisat2 command for single-end short read data.") 
        sys.exit(1)

    try:
        if not rnaseq_paired_sets == []:
            if file_format(rnaseq_paired_sets[0]) == "fasta":
                format_option = "-f"
            if file_format(rnaseq_paired_sets[0]) == "fastq":
                format_option = ""
            if file_format(rnaseq_paired_sets[0]) == "unknown":
                print("Error: Unknown file format for paired-end short reads. Please provide a FASTA or FASTQ file.")
                sys.exit(1)
            string_with_first = ",".join(rnaseq_paired_sets[0::2])
            string_with_second = ",".join(rnaseq_paired_sets[1::2])
            output_2 = "alignment_paired_rnaseq.sam"  
            hisat2_paired_command = [
                    path, 
                    format_option,                       
                    "-x", "genome",            
                    "-1", string_with_first,
                    "-2", string_with_second,   
                    "--dta",
                    "-p", str(threads),
                    "-S", output_2,                  
                ]
            print("Mapping of paired-end short reads to the genome...")
            
            result = subprocess.run(hisat2_paired_command, capture_output=True)
            
            if result.returncode == 0:
                print("Mapping of paired-end short reads completed successfully")
                alignments_list.append(output_2)
            else:
                print("Error during mapping of paired-end short reads:")
                print(result.stdout)
                print(result.stderr)
                sys.exit(1)

    except Exception:
        print("Could not run hisat2 command for paired-end short read data.")
        sys.exit(1)

    return alignments_list 

''' Mapping long reads to the genome using Minimap2 '''
def mapping_long(genome, isoseq_sets):
    try :
        output_sam = "alignment_isoseq.sam"
        path = minimap + "/minimap2" 
        command = [path,  "-ax", "splice:hq", "-uf", genome, "-t", str(threads)] + isoseq_sets + ["-o", output_sam]

        print("Mapping long reads to the genome...")

        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Mapping of long reads completed successfully.")
            return output_sam
        else:
            print("Error during mapping of long reads:")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)

    except Exception:
        print("Could not run minimap2 command.")
        sys.exit(1)

''' Converting files in SAM format into files in BAM format using SAMtools'''
def sam_to_bam(sam_file_list):
    try:  
        path = samtools + "/samtools"  
        for samfile in sam_file_list:
            output_bam = file_name(samfile) + ".bam"
            command = [
                path,
                "sort",
                "-@",
                str(threads), 
                samfile,
                "-o",
                output_bam
            ]

        print("Converting " + samfile +" to " + output_bam + "...")
        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Conversion from SAM to BAM file completed successfully.")

        else:
            print("Error during conversion:")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)

    except Exception:
        print("Could not run SAMtools command.")
        sys.exit(1)

''' Merging two BAM files into one using SAMtools'''
def merge_bam_files(bamfile_1, bamfile_2): 
    try:
        output_bam = "alignment_merged_rnaseq.bam"
        path = samtools + "/samtools"
        command = [
            path,
            "merge",
            "-f", 
            "-o",
            output_bam,
            bamfile_1,
            bamfile_2
        ]
        print("Merging BAM files " + bamfile_1 + " and " + bamfile_2 + "...")
        result = subprocess.run(command, capture_output=True) 

        if result.returncode == 0:
            print("Merged BAM files into a single BAM file successfully.")
            return output_bam
        else:
            print("Error during merging BAM files:")
            print(result.stdout)
            print(result.stderr) 
            sys.exit(1)

    except Exception:
        print("Could not run SAMtools command.")
        sys.exit(1)

''' Assembling the mapped reads using StringTie '''
def assembling(alignment_rnaseq, alignment_isoseq):
    try:
        print("Assembling the reads...")
        output_gtf = "transcripts.gtf"
        path = stringtie + "/stringtie"

        if args.mixed:
            command_mixed = [  
                path,
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
                print("Assembled short and long reads successfully.")

            else:
                print("Error during read assembly:")
                print(result.stdout)
                print(result.stderr) 
                sys.exit(1)

        if args.rnaseq:
            command_rnaseq = [
                path,
                "-p",
                str(threads),
                "-o",
                output_gtf,
                alignment_rnaseq
            ]
            result = subprocess.run(command_rnaseq, capture_output=True)

            if result.returncode == 0:
                print("Assembled short reads successfully.")
            else:
                print("Error during read assembly:")
                print(result.stdout)
                print(result.stderr)
                sys.exit(1) 
            
        if args.isoseq:
            command_isoseq = [
                path,
                "-L",
                "-p",
                str(threads),
                "-o",
                output_gtf,
                alignment_isoseq
            ]

            result = subprocess.run(command_isoseq, capture_output=True)

            if result.returncode == 0:
                print("Assembled long reads successfully")
            else:
                print("Error during read assembly:")
                print(result.stdout)
                print(result.stderr) 
                sys.exit(1)

    except Exception:
        print("Could not run stringtie command.")
        sys.exit(1)

''' Identifying ORFs within the transcripts predicted by StringTie using TransDecoder '''
def orfsearching(genome_fa, transcripts_gtf):
    try:   
        path = transdecoder2 + "/gtf_genome_to_cdna_fasta.pl"
        fasta_file_command = [
            path,
            transcripts_gtf,
            genome_fa
        ]
        print("Creating a FASTA file using the genome and the transcripts.gtf files...")

        with open("transcripts.fasta", "w") as output:
            result = subprocess.run(fasta_file_command, stdout=output, stderr=subprocess.PIPE)

        if result.returncode == 0:
            print("Created transcripts.fasta successfully.")
        else:
            print("Error during creation of transcripts.fasta file.")
            print(result.stdout)
            print(result.stderr)  
            sys.exit(1)  
    
    except Exception: 
        print("Could not run TransDecoder module gtf_genome_to_cdna_fasta.pl.")
        sys.exit(1)
    
    try:
        path = transdecoder + "/TransDecoder.LongOrfs"
        longORF_command = [
            path,
            "-t",
            "transcripts.fasta"
        ]

        print("Extracting the long open reading frames...")

        result= subprocess.run(longORF_command, capture_output=True)
        if result.returncode == 0:
            print("Extracted ORFs successfully.")
        else:
            print("Error during searching for long ORFs.")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)

    except Exception:
        print("Could not run TransDecoder.LongOrfs module.")
        sys.exit(1)

    try:
        path = transdecoder + "/TransDecoder.Predict"
        predict_command = [
            path,
            "-t",
            "transcripts.fasta"
        ]

        print("Predicting likely coding regions...")

        result= subprocess.run(predict_command, capture_output=True)
        if result.returncode == 0:
            print("Predicted likely coding regions successfully.")
        else:
            print("Error during prediction of likely coding regions.")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)

    except Exception:
        print("Could not run TransDecoder.Predict module.") 
        sys.exit(1)

def convert_gtf_to_gff3(transcripts_gtf, output_name):
    try:
        path = transdecoder2 + "/gtf_to_alignment_gff3.pl"
        command = [
            path,
            transcripts_gtf
        ]
        print("Converting GTF file " + transcripts_gtf + " into GFF3 format...")
    
        with open(output_name, "w") as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE)

        if result.returncode == 0:
            print("Converted GTF to GFF3 file successfully.")

        else:
            print("Error during converting GTF to GFF3 file.")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)
    
    except Exception:
        print("Could not run TransDecoder module gtf_to_alignment_gff3.pl.") 
        sys.exit(1)

def protein_aligning(genome, protein, alignment_scoring):
    try: 
        path = miniprot + "/miniprot"
        command = [
            path,
            "-t",      
            str(threads),
            "--genome",
            genome,
            "--protein",
            protein,
            "--aln"
        ]
        print("Aligning proteins to the genome...")
        with open("miniprot.aln", "w") as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE)
        if result.returncode == 0:
            print("Proteins aligned to the genome successfully.")
        else:
            print("Error during protein alignment with miniprot")
            print(result.stdout.decode())
            print(result.stderr.decode())
            sys.exit(1)

    except Exception:
        print("Could not run miniprot command.")
        sys.exit(1)

    try: 
        path = miniprot_boundary_scorer + "/miniprot_boundary_scorer"
        if alignment_scoring == None:
            alignment_scoring = miniprot_boundary_scorer + "/blosum62.csv"
        command = f"{path} -o miniprot_parsed.gff -s {alignment_scoring} < miniprot.aln"
        print("Scoring the protein to genome alignment...")
        
        result = subprocess.run(command, shell=True, capture_output=True)

        if result.returncode == 0:
            print("Alignment scored successfully with miniprot_boundary_scorer.")
        else:
            print("Error during scoring the alignment with miniport_boundary_scorer!")
            print(result.stdout.decode())
            print(result.stderr.decode())
            sys.exit(1)
    
    except Exception:
        print("Could not run miniprot_boundary_scorer command.")
        sys.exit(1)

''' Creating a file with the originally incomplete CDS, shortened to begin at the first start codon. '''
def shorten_incomplete_Orfs(transdecoder_pep):
    #File "shortened_candidates.pep" is used to store the shortened, originally incomplete (and now complete) ORFs
    with open("shortened_candidates.pep", "w") as output:
        #Parse the input file given by TransDecoder with all complete and incomplete ORFs
        for record in SeqIO.parse(transdecoder_pep, "fasta"):
            #ORFs with type 5prime_partial or internal are lacking a start codon
            if "type:5prime_partial" in record.description or "type:internal" in record.description:
                #Find the position of the first methionine (M) within the ORF
                m_position = record.seq.find("M")
                #If no M is found within the ORF, the position is -1
                if m_position != -1:
                    #Shorten the ORF to begin at the first M
                    record.seq = record.seq[m_position:]
                    #Find the coordinates of the old start and stop positions within the description of the incomplete ORF
                    description = record.description.split(" ")
                    coords = re.search(r":(\d+)-(\d+)\([\+\-]\)", description[7])
                    #Calculate new ORF length with the new start position. Use -1 so that the Stop * is not taken into account
                    new_length = len(record.seq) - 1 
                    if coords:
                        old_start = int(coords.group(1)) 
                        #The coordinates within the description are nucleotide coordinates, so the new start position 
                        # is calculated by multiplying the amino acid position by 3 and add that to the old position. 
                        #The m_position is 0-based, so the sum is not subtracted by 1.
                        new_start = old_start + m_position*3 
                        stop = int(coords.group(2))
                        #Replace the old description for the incomplete ORF by the new description for the complete ORF.
                        description = re.sub(f"{old_start}-{stop}", f"{new_start}-{stop}", record.description)
                        description = re.sub(r"len:\d+", f"len:{new_length}", description)
                        record.description = description
                        SeqIO.write(record, output, "fasta")
    print("Shortened ORFs that are 5' partial until the first start codon within the ORF sequence successfully.")

''' Creating a protein database from the protein input file using DIAMOND. '''
def make_diamond_db(protein_file):
    try:
        path = diamond + "/diamond"
        command = [
            path,
            "makedb",
            "--in",
            protein_file,
            "-d",
            "protein_db"
        ]
        print("Creating a protein database...")
        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Protein database created successfully.")
        else:
            print("Error during creating protein database with DIAMOND.")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)
    
    except Exception:
        print("Could not run DIAMOND makedb command.")
        sys.exit(1)

''' Searching for proteins that align with the CDS predictions made by TransDecoder using DIAMOND '''
def validating_ORFs(transdecoder_pep, output_tsv):
    try:
        path = diamond + "/diamond"
        command = [
            path,
            "blastp",
            "-d",
            "protein_db",
            "-q",
            transdecoder_pep,
            "-o",
            output_tsv
        ]
        print("Searching for matches in protein database for " + file_name(transdecoder_pep) + "...")
        result = subprocess.run(command, capture_output=True)

        if result.returncode == 0:
            print("Blastp search completed successfully and ", output_tsv, " was created.")
        else:
            print("Error during DIAMOND blastp search.")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)
    
    except Exception:
        print("Could not run DIAMOND blastp command.")
        sys.exit(1)

''' Classifying originally incomplete predictions as complete or incomplete based on the DIAMOND database search '''
def get_cds_classification(normal_tsv, shortened_tsv):
    print("Looking at the incomplete ORF predictions...")
    #Parse the DIAMOND results into a dataframe with headers equivalent to the DIAMOND output using pandas.
    #normal_tsv stored the path to the DIAMOND results for the original TransDecoder output with complete and incomplete candidates and
    #shortened_tsv stores the path to the DIAMOND results for the originally incomplete candidates shortened until the first M.
    header_list = ["cdsID", "proteinID", "percIdentMatches", "alignLength", "mismatches", "gapOpenings", "queryStart", "queryEnd", "targetStart", "targetEnd", "eValue", "bitScore"]
    df_shortened = pd.read_csv(shortened_tsv, delimiter='\t', header=None, names=header_list)
    df_normal = pd.read_csv(normal_tsv, delimiter='\t', header=None, names=header_list)
    #Merge the two dataframes on the columns "cdsID" and "proteinID" to compare the results of the two DIAMOND searches.
    merged_df = pd.merge(df_shortened, df_normal, on=["cdsID", "proteinID"], suffixes=('_short', '_normal'))
    merged_df = merged_df.drop(columns=["alignLength_short", "alignLength_normal", "mismatches_short", "mismatches_normal", "gapOpenings_short", "gapOpenings_normal", "queryEnd_short", "queryEnd_normal", "targetEnd_short", "targetEnd_normal", "eValue_short", "eValue_normal"])
    #Add a column to store the support score for each candidate.
    merged_df["supportScore"] = None
    
    #Calculate the support score for each cds-candidate.
    for i, cds in merged_df.iterrows():
        q_incomplete_start = cds["queryStart_normal"]
        t_incomplete_start = cds["targetStart_normal"]
        t_complete_start = cds["targetStart_short"] 
        aai_incomplete = cds["percIdentMatches_normal"] 
        aai_complete = cds["percIdentMatches_short"] 

        #If the AAI is 0, set it to a small value to avoid division by zero.
        if aai_complete == 0:
            aai_complete = 0.0001
        match_log = math.log(aai_incomplete/aai_complete)

        support_score = (t_complete_start - t_incomplete_start) - (q_incomplete_start - 1) + match_log**1000
        merged_df.at[i, "supportScore"] = support_score
    
    #Add a column to store the maximum bit score for each candidate, either from the shortened or normal candidate.
    merged_df["bitScore_max"] = merged_df[["bitScore_short", "bitScore_normal"]].max(axis=1) 
    #Group the dataframe by the cdsID and select the 25 best alignments based on the bitscore.
    grouped = merged_df.groupby("cdsID")
    #Create a dictionary to store the classification for each candidate.
    classifications = {}
    for groupname, groupdata in grouped:
        sorted_group = groupdata.sort_values(by="bitScore_max", ascending=False)
        incomplete = False
        #Check if the 25 best alignments contain any support score greater than 0. 
        # If so the candidate remains incomplete, otherwise it is reclassified as complete.
        for i, cds in sorted_group.head(25).iterrows():
            if cds["supportScore"] > 0:
                incomplete = True
                break
        if incomplete:
            classifications[cds["cdsID"]] = "incomplete" 
        else:
            classifications[cds["cdsID"]] = "complete"  

    print("Compared if the shortened or the original CDS has better protein support for each incomplete CDS candidate successfully.")
    return classifications
   
''' Creating revised file with translated ORFs, containing the originally and newly complete ORFs and the incomplete ORFs. '''
def get_optimized_pep_file(normal_pep, shortened_pep, classifications):
    with open("revised_candidates.pep", "w") as output:
        #Parse shortened ORFs into a dictonary 
        shortened_pep_dict = {record.id: (record.seq, record.description) for record in SeqIO.parse(shortened_pep, "fasta")}
        #Go through the original file with all ORFs 
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

            #ID can be found in the classification dictionary if it originally lacked a start codon and if any protein aligned with the ORF.
            if classification == "5prime_partial" or classification == "internal":
                if id in classifications: 
                    #If the classification is incomplete no change is made and the record from the original PEP file is written to the output file.
                    if classifications[id] == "incomplete":
                        SeqIO.write(record, output, "fasta")  
                    #Else the shortened version is written to the output file with the type replaced in the description.
                    else:
                        seq = shortened_pep_dict[id][0]
                        record.seq = seq
                        description = shortened_pep_dict[id][1]
                        if "type:5prime_partial" in description:
                            description = description.replace("type:5prime_partial", "type:complete")
                        else:
                            description = description.replace("type:internal", "type:3prime_partial")
                        record.description = description   
                        SeqIO.write(record, output, "fasta")
                #If no protein aligns with the candidate, the original record is written to the output file.
                else:
                    SeqIO.write(record, output, "fasta") 
            #If the ORF is originally complete or 3' partial, the record is written to the output file.
            else:
                SeqIO.write(record, output, "fasta")
    print("Created a revised CDS file using the classifications as complete or incomplete CDS successfully.")

''' Creating a GFF3 file from the TransDecoder PEP and the StringTie GTF files with transcript coordinates. '''
'''In training mode the GFF3 file only contains complete gene structures, in hints mode all predictions are included.'''
def from_pep_file_to_gff3(orf_pep, transcript_gtf, output_name):
    transcripts = {}
    transcript_length = 0
    with open(transcript_gtf, "r") as transcript_file:
        #Calculate the total length of each transcript by summing up the exon lengths and store it in a dictonary.
        for line in transcript_file:
            if line.startswith("#"):
                continue
            else:
                part = line.strip().split('\t')
                if "exon" == part[2]:
                    start = int(part[3])
                    stop = int(part[4])
                    length = stop - start + 1 
                    #Store the entire length of a transcript by adding up the exon lengths
                    transcript_length += length
                    transcript_id = re.search(r'transcript_id "([^"]+)"', part[8]).group(1)
                    transcripts[transcript_id] = transcript_length
                else:
                    transcript_length = 0
    #Fill the output GFF3 file with the records for each gene: mRNA, exon, CDS, 5' UTR and 3' UTR (if present).
    with open(output_name, "w") as output:
        for record in SeqIO.parse(orf_pep, "fasta"):
            id_transdecoder = record.id
            id_stringtie = id_transdecoder.split(".p")
            id_stringtie = id_stringtie[0]
            tool1 = "StringTie"
            tool2 = "TransDecoder"
            description = record.description
            coords = re.search(r":(\d+)-(\d+)\([\+\-]\)", description)
            orf_start = int(coords.group(1))
            orf_stop = int(coords.group(2))
            transcript_length = transcripts[id_stringtie]
            strand = re.search(r"\((\+|-)\)", description) 
            strand = strand.group(1) 
            description_parts = description.split()
            gene_id = description_parts[1]
            gene_name = description_parts[2] + " "+ description_parts[3] +" "+ description_parts[4]  
            output.write(f"{id_stringtie}\t{tool1}\tgene\t1\t{transcript_length}\t.\t{strand}\t.\tID={gene_id};Name=\"{gene_name}\"\n")
            output.write(f"{id_stringtie}\t{tool1}\tmRNA\t1\t{transcript_length}\t.\t{strand}\t.\tID={id_transdecoder};Parent={gene_id};Name=\"{gene_name}\"\n")
            if strand == "+" and orf_start > 1:
                output.write(f"{id_stringtie}\t{tool2}\tfive_prime_UTR\t1\t{orf_start-1}\t.\t{strand}\t.\tID={id_transdecoder}.utr5p1;Parent={id_transdecoder}\n")
            if strand == "-" and orf_stop < transcript_length:
                output.write(f"{id_stringtie}\t{tool2}\tfive_prime_UTR\t{orf_stop+1}\t{transcript_length}\t.\t{strand}\t.\tID={id_transdecoder}.utr5p1;Parent={id_transdecoder}\n")
            output.write(f"{id_stringtie}\t{tool1}\texon\t1\t{transcript_length}\t.\t{strand}\t.\tID={id_transdecoder}.exon1;Parent={id_transdecoder}\n")
            output.write(f"{id_stringtie}\t{tool2}\tCDS\t{orf_start}\t{orf_stop}\t.\t{strand}\t0\tID=cds.{id_transdecoder};Parent={id_transdecoder}\n")
            if strand == "-" and orf_start > 1:
                output.write(f"{id_stringtie}\t{tool2}\tthree_prime_UTR\t1\t{orf_start-1}\t.\t{strand}\t.\tID={id_transdecoder}.utr3p1;Parent={id_transdecoder}\n")
            if strand == "+" and orf_stop < transcript_length:
                output.write(f"{id_stringtie}\t{tool2}\tthree_prime_UTR\t{orf_stop+1}\t{transcript_length}\t.\t{strand}\t.\tID={id_transdecoder}.utr3p1;Parent={id_transdecoder}\n")
            output.write("\n")
    print("Created a GFF3 file from the " + file_name(orf_pep) + ".pep file successfully.")

''' Creating FASTA file from the records stored in dictionary. '''
def from_dict_to_pep_file(input_dict, output_name):
    with open(output_name, "w") as output:
        for entry in input_dict:
            #Create a SeqRecord using the Seq and description from the dictionary.
            record = SeqRecord(
                input_dict[entry][1],  
                id=entry,                 
                description=input_dict[entry][0]                                                           
            )
            #Write the record to the output file in FASTA format.
            SeqIO.write(record, output, "fasta")

''' Creating a file in BED format from the miniprot gene predictions that can be used to identify conflicting TransDecoder predictions. ''' 
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

''' Creating a file in BED format from the CDS candidates that can be used to find conflicts with the miniprot predictions. '''
def preparing_candidates_for_conflict_comparison(candidates_gff3, transcripts_gtf):
    #Find chromosome number for each gene in the transcripts.gtf file from StringTie and store it in a dictionary 
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

    #Create BED file
    with open("candidates.bed", "w") as output, open(candidates_gff3, "r") as candidates_file:
        for line in candidates_file:
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

''' Identifying CDS candidates that overlap with miniprot gene predictions using bedtools. '''
def finding_protein_conflicts(candidates_bed, reference_bed):
    try:
        path = bedtools + "/bedtools"
        command = [
            path,
            "coverage",
            "-a",
            candidates_bed,
            "-b",
            reference_bed,
            "-s" 
        ]
        print("Identifying overlaps between the CDS candidates and the protein-based gene predictions...")

        with open("conflicts.bed", "w") as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE)

        if result.returncode == 0:
            print("Coverage file that contains overlapping loci created successfully.")

        else:
            print("Error during creating coverage file.")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)

    except Exception:
        print("Could not run bedtools command.")
        sys.exit(1)

''' Searching for a stop codon in the 5' UTR of the CDS candidates. '''
def finding_stop_in_utr(transcripts_fasta, intrinsic_candidates_genome):
    print("Searching for stop codons in the 5' UTR of the candidates...")
    #Dictionary to store the candidates with a stop codon in the 5' UTR (True).
    stop_in_utr_dict = {}   
    #Dictionary to store the nucleotide sequence of each candidate from the transcripts.fasta  
    transcripts_dict = {} 
    #Stop codons that are searched for in the 5' UTR.                  
    stop_codons = ["TAA", "TAG", "TGA"]
    for record in SeqIO.parse(transcripts_fasta, "fasta"):
        transcripts_dict[record.id] = [record.description, record.seq]

    #Go through the candidates (provided with genome coords) 
    with open(intrinsic_candidates_genome, "r") as candidates_file:
        for line in candidates_file:
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

            #If the feature is a 5' UTR, the sequence is extracted and searched for a stop codon.
            if feature_type == "five_prime_UTR":
                utr_length = end - start + 1
                utr_sequence = transcripts_dict[stringtie_id][1][utr_length:]
                for codon in stop_codons:                       
                    stop_codon_position = utr_sequence.find(codon)
                    #If a stop codon position is found, the dictionary stores "True" for the transcript 
                    # and no other stop codon is searched for within this transcript
                    if stop_codon_position != -1:
                        stop_in_utr_dict[stringtie_id] = True
                        break
    print("Identified the stop codons in the 5' UTR of the candidates successfully.")
    return stop_in_utr_dict 

''' Selecting a set of high-confidence genes based on protein evidence. '''
def getting_hc_supported_by_proteins(diamond_tsv, transdecoder_pep, protein_file):
    print("Selecting high-confidence genes supported by protein evidence...")
    #Dictionary to store the length of each target protein sequence from the protein input file (that functions as target database)
    t_length_dict = {}
    for record in SeqIO.parse(protein_file, "fasta"):
        t_length_dict[record.id] = len(record.seq) 
    #Dictionary to store the query protein sequence and description for each CDS candidate
    q_dict = {}
    #List to store the IDs of the high-confidence genes supported by protein evidence
    already_hc_genes = []

    for record in SeqIO.parse(transdecoder_pep, "fasta"):
        q_dict[record.id] = [record.description, record.seq] 

    #Go through the DIAMOND results and check if the CDS candidates are supported by protein evidence by using the equations.
    with open(diamond_tsv, "r") as tsv, open("hc_genes.pep", "w") as output:
        for line in tsv:    
            part = line.strip().split('\t')
            id = part[0] #Query ID
            protein_id = part[1] #Target ID
            aaident = float(part[2]) #Percentage of identical amino acids 
            align_length = int(part[3]) #Length of the alignment 
            q_start = int(part[6]) #Start of the alignment in the query protein
            q_end = int(part[7]) #End of the alignment in the query protein
            t_start = int(part[8]) #Start of the alignment in the target protein
            t_end = int(part[9]) #End of the alignment in the target protein
            #Find query ID and its sequence and description in the dictionary 
            if id in q_dict:
                record.id = id 
                record.description = q_dict[id][0]
                record.seq = q_dict[id][1]
                q_length = len(record.seq) - 1  #-1 so that the stop * doesn't count
                t_length = t_length_dict[protein_id]
                
                #Only complete candidates are considered as high-confidence CDS.
                if "type:complete" in record.description: 
                    if (q_start - t_start) < 6 and (t_length - align_length) < 15 and aaident > 95:
                        SeqIO.write(record, output, "fasta")
                        gene_id = id.split(".")[0] + "." + id.split(".")[1] 
                        #If a candidate is supported by protein evidence, it is excluded from further consideration. 
                        del q_dict[id]
                        already_hc_genes.append(gene_id)

        #Only store the CDS in the dictionary that don't belong to a gene that is high-confidence with a different CDS
        q_dict = {
            id: value
            for id, value in q_dict.items()
            if id.split(".")[0] + "." + id.split(".")[1] not in already_hc_genes
        }
    print("Set of high-confidence genes was selected successfully.")
    #Return the CDS candidates that are left for further intrinsic analysis.
    return q_dict

''' Appending the set of protein supported high-confidence genes with high-confidence genes supported by intrinsic criteria. '''
def getting_hc_supported_by_intrinsic(q_dict):
    print("Selecting high-confidence genes supported by intrinsic criteria...")
    #Get a PEP file from the candidates dictionary.
    from_dict_to_pep_file(q_dict, "intrinsic_candidates.pep")
    #Select only one isoform for each gene.
    choose_one_isoform("intrinsic_candidates.pep", "intrinsic_one_isoform.pep")
    #Create a GFF3 file from the PEP file.
    from_pep_file_to_gff3("intrinsic_one_isoform.pep", "transcripts.gtf", "intrinsic_candidates.gff3")
    #Transform the transcript coordinates of the candidates into genome coordinates.
    from_transcript_to_genome("intrinsic_candidates.gff3", "transcripts.gff3", "transcripts.fasta", "intrinsic_candidates_genome.gff3")
    #Prepare a file with the candidates for conflict comparison with the miniprot predictions.
    preparing_candidates_for_conflict_comparison("intrinsic_candidates_genome.gff3", "transcripts.gtf")
    #Prepare a file with the miniprot predictions for conflict comparison with the candidates.
    preparing_miniprot_gff_for_conflict_comparison("miniprot_parsed.gff")
    #Find conflicts between the candidates and the miniprot predictions.
    finding_protein_conflicts("candidates.bed", "reference.bed")
    #Find candidates with a stop codon in the 5' UTR.
    stop_in_utr_dict= finding_stop_in_utr("transcripts.fasta", "intrinsic_candidates_genome.gff3")

    #Open bedtools output file. If a CDS candidate has no conflict with any miniprot prediction, "True" is stored in 
    # the dictionary for the transcript id and otherwise "False".
    with open("conflicts.bed", "r") as conflicts_file:
        no_conflicts_dict = {}
        for line in conflicts_file:
            part = line.split('\t')
            transcript_id = part[3]
            if int(part[6]) == 0: 
                no_conflicts_dict[transcript_id] = True
            else:
                no_conflicts_dict[transcript_id] = False

    with open("hc_genes.pep", "a") as output_hc, open("intrinsic_one_isoform.pep", "r") as candidates_file:
        already_hc_genes = []
        #Go through the candidates and decide wether they meet the criteria for high-confidence genes.
        for record in SeqIO.parse(candidates_file, "fasta"):
            transcript_id = record.id
            length_pep = int(re.search(r"len:(\d+)", record.description).group(1))
            length_cds = length_pep*3
            score = float(re.search(r"score=(-?[\d.]+)", record.description).group(1))
            if transcript_id in no_conflicts_dict: 
                condition1 = (length_cds >= 300) #Minimum 300 nucleotides in length
                condition2 = (score > 50) #Score higher than 50
                condition3 = (no_conflicts_dict[transcript_id]) #No conflicts with any miniprot prediction
                condition4 = "type:complete" in record.description #Only complete candidates are considered
                id_stringtie = transcript_id.split(".p")
                id_stringtie = id_stringtie[0]
                if id_stringtie in stop_in_utr_dict:
                    condition5 = (stop_in_utr_dict[id_stringtie]) #Must have stop codon in 5'UTR
                else:
                    condition5 = False
                #If a CDS meets the criteria it is added to the set of high-confidence CDS. 
                if condition1 and condition2 and condition3 and condition4 and condition5:
                    id_stringtie = transcript_id.split(".p")
                    id_stringtie = id_stringtie[0]
                    already_hc_genes.append(id_stringtie)
                    SeqIO.write(record, output_hc, "fasta")
        print("Appended set of high-confidence genes with genes supported by intrinsic criteria successfully.")

        #File with low-confidence genes is filled: Every gene that contains no high-confidence CDS is added
        with open("intrinsic_one_isoform.pep", "r") as candidates_file, open("lc_genes.pep", "w") as output_lc:
            for record in SeqIO.parse(candidates_file, "fasta"):
                transcript_id = record.id
                id_stringtie = transcript_id.split(".p")
                id_stringtie = id_stringtie[0]  
                if id_stringtie not in already_hc_genes:
                    SeqIO.write(record, output_lc, "fasta")
            print("All not-high-confidence genes create the set of low-confidence genes.")

''' Select only one isoform for each gene. '''
def choose_one_isoform(transdecoder_pep, output_name):
    #Dictionary to store a list of isoforms for each gene.
    isoform_dict = {}
    #Write a new file with just one isoform per a gene.
    with open(output_name, "w") as output:
        #First, all isoforms of each gene are stored in a list in the dictionary with the gene id as key. 
        for record in SeqIO.parse(transdecoder_pep, "fasta"):
            transdecoder_id = record.id
            stringtie_id = transdecoder_id.split(".p")[0]
            gene_id = stringtie_id.split(".")[1]
            description = record.description.split(" ")
            cds_coords = re.search(r":(\d+)-(\d+)\((\+|\-)\)", description[7])
            start_cds = int(cds_coords.group(1))
            stop_cds = int(cds_coords.group(2))
            if gene_id not in isoform_dict:
                isoform_dict[gene_id] = [(start_cds, stop_cds, record)]
            else:
                isoform_dict[gene_id].append((start_cds, stop_cds, record))

        #Secondly, one isoform is selected for each gene and written to the output file.
        for gene_id in isoform_dict:
            if len(isoform_dict[gene_id]) == 1:
                record = isoform_dict[gene_id][0][2]
                SeqIO.write(record, output, "fasta")
            else:
                #The longest isoform is chosen. 
                #If there are two longest, the first one is taken, because that indicates the highest expression rate. 
                longest_isoform = max(isoform_dict[gene_id], key=lambda x: x[1] - x[0]) 
                record = longest_isoform[2]
                SeqIO.write(record, output, "fasta")
    print("Selected one isoform for each gene successfully.") 

''' StringTie provides transcript and exon boundaries in genome coordinates, while TransDecoder provides CDS boundaries in transcript coordinates. 
    This function transforms the transcript coordinates into genome coordinates using a TransDecoder module. '''
def from_transcript_to_genome(cds_gff3, transcripts_gff3, transcripts_fasta, output_name):
    try:
        path = transdecoder2 + "/cdna_alignment_orf_to_genome_orf.pl"
        command = [
            path,
            cds_gff3,
            transcripts_gff3,
            transcripts_fasta
        ]
        print("Transforming transcript coordinates into genome coordinates of " + file_name(cds_gff3) + "...")
        with open(output_name, "w") as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE)

        if result.returncode == 0:
            print("Transformed transcript into genome coordinates successfully.")

        else:
            print("Error during the transformation of transcript coordinates into genome coordinates:")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)
    
    except Exception:
        print("Could not run TransDecoder module cdna_alignment_orf_to_genome_orf.pl")
        sys.exit(1)

        
''' Creating a hints file with CDS, start and stop codon and intron hints, using the high-confidence and low-confidence files with 
    genome coordinates that contains mRNA, exon, CDS and UTR features. '''
def creating_intron_hints_file(hints_hc_genome_gff3, hints_lc_genome_gff3, output_name):
    with open(hints_hc_genome_gff3, "r") as hc_file, open(hints_lc_genome_gff3, "r") as lc_file, open(output_name, "w") as output:
        #List to store exon and CDS features in order determine the introns and start/stop codons.
        exon_list = []
        cds_list = []
        for file in [hc_file, lc_file]:
            if file == hc_file:
                #Confidence level which is added to the attributes of the features. M indicated a high confidence and P a low confidence gene
                confidence_level = "M"
            if file == lc_file:
                confidence_level = "P"
            for line in file:
                part = line.strip().split('\t')
                if len(part) < 9: #Skip empty lines
                    continue
                feature = part[2]
                #If a new mRNA is listed, the intron and start/stop features of the last mRNA are added to the file.
                if feature == "mRNA":
                    if len(exon_list) > 1:
                        strand = exon_list[0][6]
                        if strand == "+":
                            for i in range(len(exon_list) - 1):
                                current_exon = exon_list[i]
                                next_exon = exon_list[i + 1]
                                intron_start = int(current_exon[4]) + 1
                                intron_stop = int(next_exon[3]) - 1
                                chromosome = current_exon[0]
                                attributes = current_exon[8]
                                updated_attributes = re.sub(r'exon(\d+)', r'intron\1', attributes)
                                updated_attributes_with_level =  "src="+ confidence_level + ";"+updated_attributes
                                output.write(f"{chromosome}\tStringTie\tintron\t{intron_start}\t{intron_stop}\t.\t{strand}\t.\t{updated_attributes_with_level}\n")
                        if strand == "-":
                            for i in range(len(exon_list) - 1, 0, -1):
                                current_exon = exon_list[i]
                                next_exon = exon_list[i - 1]
                                intron_start = int(current_exon[4]) + 1
                                intron_stop = int(next_exon[3]) - 1
                                chromosome = current_exon[0]
                                attributes = current_exon[8]
                                intron_number = "intron" + str(i)
                                updated_attributes = re.sub(r'exon(\d+)', intron_number, attributes)
                                updated_attributes_with_level =  "src="+ confidence_level + ";"+updated_attributes 
                                output.write(f"{chromosome}\tStringTie\tintron\t{intron_start}\t{intron_stop}\t.\t{strand}\t.\t{updated_attributes_with_level}\n")
                    if len(cds_list) > 0:
                        first_cds = cds_list[0]
                        start_codon_start = first_cds[3]
                        start_codon_stop = int(first_cds[3]) + 2
                        last_cds = cds_list[-1]
                        stop_codon_start = int(last_cds[4]) - 2
                        stop_codon_stop = last_cds[4]
                        chromosome = first_cds[0]
                        strand = first_cds[6]
                        attributes_first = first_cds[8]
                        attributes_last = last_cds[8]
                        attributes_first_with_level =  "src="+ confidence_level + ";" + attributes_first
                        attributes_last_with_level =  "src="+ confidence_level + ";" + attributes_last
                        output.write(f"{chromosome}\tTransDecoder\tstart\t{start_codon_start}\t{start_codon_stop}\t.\t{strand}\t.\t{attributes_first_with_level}\n")
                        output.write(f"{chromosome}\tTransDecoder\tstop\t{stop_codon_start}\t{stop_codon_stop}\t.\t{strand}\t.\t{attributes_last_with_level}\n")
                    exon_list.clear()
                    cds_list.clear()
                else:
                    if feature == "exon":
                        #For each mRNA an exon list temporarily stores the exon coordinates until they are written to the file
                        exon_list.append(part)
                    
                    if feature == "CDS":
                        cds_list.append(part)
                        updated_attributes = "src=" + confidence_level + ";" + part[8]
                        updated_line = re.sub(part[8], updated_attributes, line)
                        #Add CDS to the file
                        output.write(updated_line)
    print("Created a hints file with intron hints successfully.")

''' Opening the config file with the input data. '''
def load_config(config_file):
    with open(config_file, "r") as config_file:
        input_files = yaml.safe_load(config_file)
        return input_files

''' CONFIGURATIONS '''
print("*************************************************************************************")
print("                                     PreGalba                                        ")
print("*************************************************************************************")
parser = argparse.ArgumentParser(description='Genome annotation with protein and transcriptomic evidence auch as RNA-Seq and Iso-Seq reads')  
parser.add_argument('-t', '--threads', default=4, help='Number of threads (default=4)', required=False)
parser.add_argument('-y', '--config', help='Config file path', metavar='<config.yaml>', required=False) 
parser.add_argument('-g', '--genome', help='Genome file path', metavar='<genome.fasta>', required=False)
parser.add_argument('-p', '--proteins', help='Protein file path', metavar='<proteins.fasta>', required=False)
parser.add_argument('-sp', '--short_paired', help='Comma separated paired-end short read file paths', metavar='<set1.fasta,set2.fasta>', required=False)
parser.add_argument('-ss', '--short_single', help='Comma separated single-end short read file paths', metavar='<set1.fasta,set2.fasta>', required=False)
parser.add_argument('-l', '--long', help='Comma separated long read file paths', metavar='<set1.fasta,set2.fasta>', required=False)
parser.add_argument('-m', '--scoring_matrix', help='Alignment scoring matrix for amino acids', metavar='<blosum62.csv>', required=False)

parser.add_argument('--isoseq', action='store_true', help='Use this option if you want to process isoseq data only')
parser.add_argument('--rnaseq', action='store_true', help='Use this option if you want to process rnaseq data only')
parser.add_argument('--mixed', action='store_true', help='Use this option if you want to process both rnaseq and isoseq data')
parser.add_argument('--projname', help='Name the output folder', required=False)
parser.add_argument('--output_path',  help='Specify the path you want the output folder to be created in if it should differ from the current working directory', required=False)

parser.add_argument('--HISAT', help='Provide the path to the HISAT2 directory', required=False)
parser.add_argument('--MINIMAP', help='Provide the path to the Minimap2 directory', required=False)
parser.add_argument('--STRINGTIE', help='Provide the path to the StringTie2 directory', required=False)
parser.add_argument('--SAMTOOLS', help='Provide the path to the samtools directory', required=False)
parser.add_argument('--TRANSDECODER', help='Provide the path to the TransDecoder directory', required=False)
parser.add_argument('--DIAMOND', help='Provide the path to the DIAMOND directory', required=False)
parser.add_argument('--BEDTOOLS', help='Provide the path to the bedtools2 directory', required=False)
parser.add_argument('--MINIPROT', help='Provide the path to the miniprot directory', required=False)
parser.add_argument('--MINIPROT_BOUNDARY_SCORER', help='Provide the path to the miniprot-boundary-scorer directory', required=False)

args = parser.parse_args()

genome_file = None
rnaseq_paired_sets = []
rnaseq_single_sets = []
isoseq_sets = []
protein_file = None
scoring_matrix = None
projname = "pregalba"

#Option 1: Providing the data in a config file in YAML format
if args.config:
    input_files = load_config(args.config)
    species_name = input_files.get("species", None)
    if species_name != None:
        projname = species_name + "_pregalba"

    genome_file = input_files.get("genome", None)
    rnaseq_paired_sets = input_files.get("rnaseq_paired_sets", []) 
    rnaseq_single_sets = input_files.get("rnaseq_single_sets", [])
    isoseq_sets = input_files.get("isoseq_sets", [])
    protein_file = input_files.get("protein", None) 
    scoring_matrix = input_files.get("scoring_matrix", None)

    hisat = input_files.get("hisat", None) 
    minimap = input_files.get("minimap", None) 
    stringtie = input_files.get("stringtie", None)
    samtools = input_files.get("samtools", None)
    transdecoder = input_files.get("transdecoder", None)
    diamond = input_files.get("diamond", None)
    bedtools = input_files.get("bedtools", None)
    miniprot = input_files.get("miniprot", None)
    miniprot_boundary_scorer = input_files.get("miniprot_boundary_scorer", None) 

else:
    scoring_matrix = None
    projname = "pregalba"

#Option 2: Providing the data with the program call
#If both is given, the program will use the data provided with the program call

if args.genome:
    genome_file = args.genome

if args.proteins:
    protein_file = args.proteins

if args.short_paired:
    rnaseq_paired_sets = args.short_paired.split(",")

if args.short_single:
    rnaseq_single_sets = args.short_single.split(",")

if args.long:
    isoseq_sets = args.long.split(",")

if args.scoring_matrix:
    scoring_matrix = args.scoring_matrix

threads = args.threads

if args.projname:
    projname = args.projname
    print("P R O J E C T: " + projname)

if args.output_path:
    output_path = args.output_path
    print("Creating output directory in " + output_path + ".")
else:
    output_path = os.getcwd()
    print("Creating output directory in the current working directory.")

os.makedirs(output_path + "/" + projname, exist_ok=True)
os.chdir(output_path + "/" + projname)

#Exit and/or print an error message if required data is missing, error-prone or incompatible with the selected mode.
if genome_file == None:
    print("Error: No genome file path found. Please provide a genome file.")
    sys.exit(1) 
if file_format(genome_file) != "fasta" :
    print("Error: Incompatible file format for genome file. Please provide a file in FASTA format.")
    sys.exit(1)
if protein_file == None:
    print("Error: No protein file path found. Please provide a protein file.")
    sys.exit(1)
if file_format(protein_file) != "fasta":
    print("Error: Incompatible file format for protein file. Please provide a file in FASTA format.")
    sys.exit(1)
if rnaseq_paired_sets == [] and rnaseq_single_sets == [] and isoseq_sets == []:
    print("Error: No transcriptomic data found in config file. Please provide at least one set of RNA-Seq or Iso-Seq data.")
    sys.exit(1)
if not args.rnaseq and not args.isoseq and not args.mixed:
    if (rnaseq_paired_sets != [] or rnaseq_single_sets != []) and isoseq_sets == []:
        args.rnaseq = True
    if (rnaseq_paired_sets == [] and rnaseq_single_sets == []) and isoseq_sets != []:
        args.isoseq = True
    if (rnaseq_paired_sets != [] or rnaseq_single_sets != []) and isoseq_sets != []:
        args.mixed = True
    print("You did not specify which data you want to process. The mode is set based on the data provided.")
if (args.rnaseq and args.isoseq) or (args.rnaseq and args.mixed) or (args.isoseq and args.mixed):
    print("You selected multiple modes. The mode was set to mixed.")
if rnaseq_paired_sets == [] and rnaseq_single_sets == [] and args.rnaseq:
    print("Error: No RNA-Seq data found in config file. Please provide at least one set of RNA-Seq data.")
    sys.exit(1)
if isoseq_sets == [] and args.isoseq:
    print("Error: No Iso-Seq data found in config file. Please provide at least one set of Iso-Seq data.")
    sys.exit(1)
if (rnaseq_paired_sets == [] and rnaseq_single_sets == []) or (isoseq_sets == []) and args.mixed:
    print("Error: You chose the mixed option. Please provide both RNA-Seq and Iso-Seq data.")
    sys.exit(1)

process_rnaseq = args.rnaseq or args.mixed
process_isoseq = args.isoseq or args.mixed

if args.rnaseq:
    mode = "rnaseq"
if args.isoseq:
    mode = "isoseq"
if args.mixed:
    mode = "mixed"

if args.HISAT:
    hisat = args.HISAT 

if args.HISAT or (args.config and hisat != None):
    find_hisat = hisat + "/hisat2"
    if not (os.access(find_hisat, os.X_OK)):
        print("Error: HISAT2 is not found or not executable at given path!")
        exit(1)
    else:
        print("HISAT2 path was set to: " + hisat)
else:
    if shutil.which("hisat2") is not None:
        hisat =  os.path.dirname(shutil.which("hisat2"))

    else:
        print("Error: hisat2 file wasn't found.")
        print("Please provide the path to the StringTie2 directory in your config file or use option --HISAT2.")
        print("Also you can use the docker file available.")
        exit(1)

if args.MINIMAP:
    minimap = args.MINIMAP 

if args.MINIMAP or (args.config and minimap != None):
    find_minimap = minimap + "/minimap2"
    if not (os.access(find_minimap, os.X_OK)):
        print("Error: Minimap2 is not found or not executable at given path!")
        exit(1)
    else:
        print("Minimap2 path was set to: " + minimap)
else:
    if shutil.which("minimap2") is not None:
        minimap =  os.path.dirname(shutil.which("minimap2"))

    else:
        print("Error: minimap2 file wasn't found.")
        print("Please provide the path to the Minimap2 directory in your config file or use option --Minimap2.")
        print("Also you can use the docker file available.")
        exit(1)

if args.STRINGTIE:
    stringtie = args.STRINGTIE 

if args.STRINGTIE or (args.config and stringtie != None):
    find_stringtie = stringtie + "/stringtie"
    if not (os.access(find_stringtie, os.X_OK)):
        print("Error: StringTie2 is not found or not executable at given path!")
        exit(1)
    else:
        print("StringTie2 path was set to: " + stringtie)
else:
    if shutil.which("stringtie") is not None:
        stringtie =  os.path.dirname(shutil.which("stringtie"))

    else:
        print("Error: stringtie file wasn't found.")
        print("Please provide the path to the StringTie2 directory in your config file or use option --STRINGTIE.")
        print("Also you can use the docker file available.")
        exit(1)

if args.SAMTOOLS:
    samtools = args.SAMTOOLS

if args.SAMTOOLS or (args.config and samtools != None):
    find_samtools = samtools + "/samtools"
    if not (os.access(find_samtools, os.X_OK)):
        print("Error: samtools is not found or not executable at given path!")
        exit(1)
    else:
        print("Samtools path was set to: " + samtools)
else:
    if shutil.which("samtools") is not None:
        samtools =  os.path.dirname(shutil.which("samtools"))

    else:
        print("Error: samtools file wasn't found.")
        print("Please provide the path to the samtools directory in your config file or use option --SAMTOOLS.")
        print("Also you can use the docker file available.")
        exit(1)

if args.TRANSDECODER:
    transdecoder = args.TRANSDECODER

if args.TRANSDECODER or (args.config and transdecoder != None):
    find_transdecoder1 = transdecoder + "/TransDecoder.LongOrfs"
    if not (os.access(find_transdecoder1, os.X_OK)):
        print("Error: Neccessary TransDecoder modules are not found or not executable at given path!")
        exit(1)
    else:
        print("TransDecoder path was set to: " + transdecoder)
        transdecoder2 = transdecoder + "/util"
else:
    if shutil.which("TransDecoder.LongOrfs") is not None:
        transdecoder =  os.path.dirname(shutil.which("TransDecoder.LongOrfs"))
        transdecoder2 =  os.path.dirname(shutil.which("gtf_genome_to_cdna_fasta.pl"))

    else:
        print("Error: TransDecoder modules weren't found.")
        print("Please provide the path to the TransDecoder directory in your config file or use option --TRANSDECODER.")
        print("Also you can use the docker file available.")
        exit(1)

if args.DIAMOND:
    diamond = args.DIAMOND 

if args.DIAMOND or (args.config and diamond != None):
    find_diamond = diamond + "/diamond"
    if not (os.access(find_diamond, os.X_OK)):
        print("Error: DIAMOND is not found or not executable at given path!")
        exit(1)
    else:
        print("DIAMOND path was set to: " + diamond)
else:
    if shutil.which("diamond") is not None:
        diamond =  os.path.dirname(shutil.which("diamond"))

    else:
        print("Error: DIAMOND file wasn't found.")
        print("Please provide the path to the DIAMOND directory in your config file or use option --DIAMOND.")
        print("Also you can use the docker file available.")
        exit(1)

if args.BEDTOOLS:
    bedtools = args.BEDTOOLS

if args.BEDTOOLS or (args.config and bedtools != None):
    find_bedtools = bedtools + "/bin/bedtools"
    if not (os.access(find_bedtools, os.X_OK)):
        print("Error: bedtools is not found or not executable at given path!")
        print("If you provided the path to the bin directory: ")
        print("Please provide the path to the bedtools directory, not to the bin directory within bedtools!")
        exit(1)
    else:
        bedtools = bedtools + "/bin"
        print("Bedtools path was set to: " + bedtools)
else:
    if shutil.which("bedtools") is not None:
        bedtools =  os.path.dirname(shutil.which("bedtools"))

    else:
        print("Error: bedtools wasn't found.")
        print("Please provide the path to the bedtools directory in your config file or use option --BEDTOOLS.")
        print("Also you can use the docker file available.")
        exit(1)

if args.MINIPROT:
    miniprot = args.MINIPROT

if args.MINIPROT or (args.config and miniprot != None):
    find_miniprot = miniprot + "/miniprot"
    if not (os.access(find_miniprot, os.X_OK)):
        print("Error: miniprot is not found or not executable at given path!")
        exit(1)
    else:
        print("Miniprot path was set to: " + miniprot)
else:
    if shutil.which("miniprot") is not None:
        miniprot =  os.path.dirname(shutil.which("miniprot"))

    else:
        print("Error: miniprot file wasn't found.")
        print("Please provide the path to the miniprot directory in your config file or use option --MINIPROT.")
        print("Also you can use the docker file available.")
        exit(1)
    
if args.MINIPROT_BOUNDARY_SCORER:
    miniprot_boundary_scorer = args.MINIPROT_BOUNDARY_SCORER

if args.MINIPROT_BOUNDARY_SCORER or (args.config and miniprot_boundary_scorer != None):
    find_miniprot_boundary_scorer = miniprot_boundary_scorer + "/miniprot_boundary_scorer"
    if not (os.access(find_miniprot_boundary_scorer, os.X_OK)):
        print("Error: miniprot-boundary-scorer is not found or not executable at given path!")
        exit(1)
    else:
        print("Miniprot-boundary-scorer path was set to: " + miniprot_boundary_scorer)
else:
    if shutil.which("miniprot_boundary_scorer") is not None:
        miniprot_boundary_scorer =  os.path.dirname(shutil.which("miniprot_boundary_scorer"))

    else:
        print("Error: miniprot_boundary_scorer file wasn't found.")
        print("Please provide the path to the miniprot_boundary_scorer directory in your config file or use option --MINIPROT_BOUNDARY_SCORER.")
        print("Also you can use the docker file available.")
        exit(1)

''' MAIN '''
print("                                                                             ")
print("Starting genome annotation in " + mode + " mode:")
print("                                                                             ")

if process_rnaseq:
    indexing(genome_file)
    alignments_list = mapping_short(rnaseq_paired_sets, rnaseq_single_sets)
    sam_to_bam(alignments_list)    
    if len(alignments_list) > 1:    
        alignment_rnaseq = file_name(merge_bam_files(alignments_list[0], alignments_list[1])) + ".bam"
    else:
        alignment_rnaseq = file_name(alignments_list[0]) + ".bam"
else:
    alignment_rnaseq = None

if process_isoseq:
    mapping_long(genome_file, isoseq_sets)
    sam_to_bam(["alignment_isoseq.sam"]) 
    alignment_isoseq = "alignment_isoseq.bam"
else:
    alignment_isoseq = None

assembling(alignment_rnaseq, alignment_isoseq)
convert_gtf_to_gff3("transcripts.gtf", "transcripts.gff3") 
orfsearching(genome_file, "transcripts.gtf")  
shorten_incomplete_Orfs("transcripts.fasta.transdecoder.pep")
make_diamond_db(protein_file)
validating_ORFs("shortened_candidates.pep", "diamond_shortened.tsv")
make_diamond_db(protein_file)
validating_ORFs("transcripts.fasta.transdecoder.pep", "diamond_normal.tsv")
classifications_dict = get_cds_classification("diamond_normal.tsv", "diamond_shortened.tsv")
get_optimized_pep_file("transcripts.fasta.transdecoder.pep", "shortened_candidates.pep", classifications_dict)
make_diamond_db(protein_file)
validating_ORFs("revised_candidates.pep", "diamond_revised.tsv")

#Getting first final output: high-confidence training.gff3
q_dict = getting_hc_supported_by_proteins("diamond_revised.tsv", "revised_candidates.pep", protein_file)
protein_aligning(genome_file, protein_file, scoring_matrix) 
getting_hc_supported_by_intrinsic(q_dict)
choose_one_isoform("hc_genes.pep", "one_chosen_isoform.pep")
from_pep_file_to_gff3("one_chosen_isoform.pep", "transcripts.gtf", "one_chosen_isoform.gff3")
from_transcript_to_genome("one_chosen_isoform.gff3","transcripts.gff3","transcripts.fasta", "training.gff3")
print("------> Training genes file training.gff3 was generated successfully!")

#Getting second final output: hints.gff3
from_pep_file_to_gff3("hc_genes.pep", "transcripts.gtf", "hc_genes.gff3")
from_pep_file_to_gff3("lc_genes.pep", "transcripts.gtf", "lc_genes.gff3")
from_transcript_to_genome("hc_genes.gff3","transcripts.gff3","transcripts.fasta", "pre_hc_hints.gff3")
from_transcript_to_genome("lc_genes.gff3","transcripts.gff3","transcripts.fasta", "pre_lc_hints.gff3")
creating_intron_hints_file("pre_hc_hints.gff3", "pre_lc_hints.gff3",  "hints.gff3")
print("------> Hints file hints.gff3 was generated successfully!")

print("                                                                                     ")
print("                                     Finished                                        ")
print("*************************************************************************************")