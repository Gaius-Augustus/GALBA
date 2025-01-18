#!/bin/bash
#SBATCH -J my_job                            # job name
#SBATCH -N 1                                  # number of nodes
#SBATCH -n 72                                  # number of MPI processes, here 1 per node
#SBATCH --partition=pinky                 # choose nodes from partition
#SBATCH -o %j.out                            # stdout file name (%j: job ID)
#SBATCH -e %j.err                             # stderr file name (%j: job ID)
#SBATCH -t 72:00:00                        # max run time (hh:mm:ss), max 72h!
#SBATCH --mail-type=end
#SBATCH --mail-user=s-amknut@uni-greifswald.de

## optional environment variables
echo "On which nodes it executes:"
echo $SLURM_JOB_NODELIST
echo " "
echo "jobname: $SLURM_JOB_NAME"

## load modules
module load singularity
#module load gcc/5.2.0   # for instance

singularity build --force galba.sif docker://amreiknuth/galba-notebook:devel 

#./scripts/main.py -g /home/s-amknut/GALBA/scripts/testGenome.fa -r /home/s-amknut/GALBA/scripts/shorttestrna1.fq
#/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/SRR12076896_1.fastq
#/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked
#./scripts/main.py -g /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked -r /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/SRR12076896_1.fastq
#singularity exec -B $PWD:$PWD galba.sif /home/s-amknut/GALBA/scripts/main.py -g /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked -r /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/SRR12076896_1.fastq
#singularity exec -B $PWD:$PWD galba.sif /home/s-amknut/GALBA/scripts/main.py -g /home/s-amknut/GALBA/entireGenome.fasta -s /home/s-amknut/GALBA/entireRna1.fastq -l /home/s-amknut/GALBA/entireLongRna.fsa
#./scripts/main.py -y arabidopsis.yaml
# ./braker.pl --species=Arabidopsis_thaliana --genome=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked --rnaseq_sets_ids=SRR12076896_1 --rnaseq_sets_dirs=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/
#WD=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana
#WD=/home/nas-hs/projs/galba-isoseq/data/Mus_musculus
#singularity exec -B $WD galba.sif /home/s-amknut/GALBA/scripts/main.py -y /home/s-amknut/GALBA/arabidopsis.yaml -t 72 --mixed   
#singularity exec -B $WD galba.sif /home/s-amknut/GALBA/scripts/main.py  -t 72 --mixed --BEDTOOLS /home/s-amknut/GALBA/tools/bedtools2 --projname pregalba_data_in_call -l /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/isoseq/m54053_170316_122436_GeneBank.fsa,/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/isoseq/m54053_170316_203815_GeneBank.fsa -sp /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/SRR12076896_1.fastq,/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/SRR12076896_2.fastq -g /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked -p /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/proteins/proteins.fa -m /home/s-amknut/GALBA/tools/blosum62_1.csv #Job 5202587
#singularity exec -B $WD galba.sif /home/s-amknut/GALBA/scripts/main.py -y /home/s-amknut/GALBA/arabidopsis.yaml -t 72 --isoseq  --BEDTOOLS /home/s-amknut/GALBA/tools/bedtools2 --projname pregalba_merge_rnaseq
#singularity exec -B $WD galba.sif /home/s-amknut/GALBA/scripts/main.py -y /home/s-amknut/GALBA/arabidopsis.yaml -t 72 --mixed --BEDTOOLS /home/s-amknut/GALBA/tools/bedtools2 --HISAT /home/s-amknut/GALBA/tools/hisat2 --MINIMAP /home/s-amknut/GALBA/tools/minimap2 --TRANSDECODER /home/s-amknut/GALBA/tools/eviann --STRINGTIE /home/s-amknut/GALBA/tools/stringtie2 --MINIPROT /home/s-amknut/GALBA/tools/miniprot  --MINIPROT_BOUNDARY_SCORER /home/s-amknut/GALBA/tools/miniprot-boundary-scorer --projname pregalba_tools_in_run
#singularity exec -B $WD galba.sif /home/s-amknut/GALBA/scripts/main.py -y /home/s-amknut/GALBA/arabidopsis.yaml -t 72 --mixed --projname pregalba_tools_in_yaml --output_path /home/s-amknut/GALBA/tools/
#singularity exec -B $WD galba.sif /home/s-amknut/GALBA/scripts/main.py -t 72 -y /home/s-amknut/GALBA/musmusculus.yaml --isoseq --BEDTOOLS /home/s-amknut/GALBA/tools/bedtools2 --projname pregalba_mus_mixed_hc
#Command wurde ausgeführt und hat braker directory produziert:
#singularity exec -B $WD /home/s-amknut/GALBA/tools/BRAKER/braker3.sif /home/s-amknut/GALBA/tools/BRAKER/scripts/braker.pl --species=Arabidopsis_thaliana --useexisting --genome=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked --rnaseq_sets_ids=SRR12076896_1,SRR12076896_2,SRR12547664_1,SRR12547664_2,SRR4010853_1,SRR4010853_2,SRR7289569_1,SRR7289569_2,SRR8714016_1,SRR8714016_2,SRR8759751_1,SRR8759751_2 --rnaseq_sets_dirs=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/ --BAMTOOLS_PATH=/home/s-amknut/GALBA/tools/bamtools/build/src/ --PROTHINT_PATH=/home/s-amknut/GALBA/tools/ProtHint/bin/ --prot_seq=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/proteins/proteins.fa 
#singularity exec -B $WD /home/s-amknut/GALBA/tools/BRAKER/braker3.sif /home/s-amknut/GALBA/tools/BRAKER/scripts/braker.pl --species=Mus_musculus_rnaseq3 --genome=/home/nas-hs/projs/galba-isoseq/data/Mus_musculus/genome/genome.fasta.masked --rnaseq_sets_ids=ERR3005082_1,ERR3005082_2,SRR10115888_1,SRR10115888_2,SRR3094250_1,SRR3094250_2,SRR5197958_1,SRR5197958_2,SRR6067921_1,SRR6067921_2,SRR9202226_1,SRR9202226_2 --rnaseq_sets_dirs=/home/nas-hs/projs/galba-isoseq/data/Mus_musculus/rnaseq/ --threads=72 --prot_seq=/home/nas-hs/projs/galba-isoseq/data/Mus_musculus/proteins/proteins.fa --bam=/home/s-amknut/GALBA/pregalba_mus_rnaseq_transdecoder/alignment_paired_rnaseq.bam --workingdir=/home/s-amknut/GALBA/braker_mus_rnaseq --skipOptimize
#singularity exec -B $WD /home/s-amknut/GALBA/tools/BRAKER/braker3.sif /home/s-amknut/GALBA/tools/BRAKER/scripts/braker.pl --species=Mus_musculus_rnaseq4 --genome=/home/nas-hs/projs/galba-isoseq/data/Mus_musculus/genome/genome.fasta.masked --threads=256 --prot_seq=/home/nas-hs/projs/galba-isoseq/data/Mus_musculus/proteins/proteins.fa --bam=/home/s-amknut/GALBA/pregalba_mus_rnaseq_transdecoder/alignment_paired_rnaseq.bam --workingdir=/home/s-amknut/GALBA/braker_mus_rnaseq 
#singularity exec -B $WD /home/s-amknut/GALBA/tools/BRAKER/braker3_lr.sif /home/s-amknut/GALBA/tools/BRAKER/scripts/braker.pl --species=Arabidopsis_thaliana_isoseq2 --genome=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked --prot_seq=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/proteins/proteins.fa --bam=/home/s-amknut/GALBA/pregalba_isoseq_transdecoder/alignment_isoseq.bam --threads=72 --workingdir=/home/s-amknut/GALBA/braker_ara_isoseq

#singularity exec -B $WD /home/s-amknut/GALBA/tools/BRAKER/braker3_lr.sif /home/s-amknut/GALBA/tools/BRAKER/scripts/braker.pl --species=Mus_Musculus_isoseq4 --genome=/home/nas-hs/projs/galba-isoseq/data/Mus_musculus/genome/genome.fasta.masked  --prot_seq=/home/nas-hs/projs/galba-isoseq/data/Mus_musculus/proteins/proteins.fa --bam=/home/s-amknut/GALBA/pregalba_mus_isoseq_transdecoder/alignment_isoseq.bam --threads=72 --workingdir=/home/s-amknut/GALBA/braker_mus_isoseq 

#RUN GENEMARK ALONE
#singularity exec -B $WD /home/s-amknut/GALBA/tools/GeneMark-ETP/bin/gmetp.pl --cores 16 --cfg /home/s-amknut/GALBA/musmusculus.yaml
