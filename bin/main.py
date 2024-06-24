#!/usr/bin/env python3

import subprocess
import sys
import os

genome_fasta = sys.argv[1]
output = sys.argv[2]
reads_file = sys.argv[3]

if not os.path.isfile(genome_fasta):
    print("Error: The file {genome_fasta} does not exist.")

elif not os.path.isfile(reads_file):
    print("Error: The file {reads_file} does not exist.")

else:
    try:
        print("Hisat2 build AMREI")
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

    try:
        print("Hisat2 AMREI")
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


        