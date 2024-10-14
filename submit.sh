#!/bin/bash
#SBATCH -J my_job                            # job name
#SBATCH -N 4                                  # number of nodes
#SBATCH -n 4                                   # number of MPI processes, here 1 per node
#SBATCH --partition=batch                  # choose nodes from partition
#SBATCH -o %j.out                            # stdout file name (%j: job ID)
#SBATCH -e %j.err                             # stderr file name (%j: job ID)
#SBATCH -t 24:00:00                        # max run time (hh:mm:ss), max 72h!
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

#./scripts/main.py -g /home/s-amknut/GALBA/scripts/testGenome.fa -r /home/s-amknut/GALBA/scripts/shorttestrna1.fq
#/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/SRR12076896_1.fastq
#/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked
#./scripts/main.py -g /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked -r /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/SRR12076896_1.fastq
#singularity exec -B $PWD:$PWD galba.sif /home/s-amknut/GALBA/scripts/main.py -g /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked -r /home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/SRR12076896_1.fastq
#singularity exec -B $PWD:$PWD galba.sif /home/s-amknut/GALBA/scripts/main.py -g /home/s-amknut/GALBA/entireGenome.fasta -s /home/s-amknut/GALBA/entireRna1.fastq -l /home/s-amknut/GALBA/entireLongRna.fsa
#./scripts/main.py -y arabidopsis.yaml
# ./braker.pl --species=Arabidopsis_thaliana --genome=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked --rnaseq_sets_ids=SRR12076896_1 --rnaseq_sets_dirs=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/
WD=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana
singularity exec -B $WD galba.sif /home/s-amknut/GALBA/scripts/main.py -t 16 -y /home/s-amknut/GALBA/arabidopsis.yaml --mixed

#Command wurde ausgeführt und hat braker directory produziert:
#singularity exec -B $WD /home/s-amknut/GALBA/tools/BRAKER/braker3.sif /home/s-amknut/GALBA/tools/BRAKER/scripts/braker.pl --species=Arabidopsis_thaliana --useexisting --genome=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/genome/genome.fasta.masked --rnaseq_sets_ids=SRR12076896_1,SRR12076896_2,SRR12547664_1,SRR12547664_2,SRR4010853_1,SRR4010853_2,SRR7289569_1,SRR7289569_2,SRR8714016_1,SRR8714016_2,SRR8759751_1,SRR8759751_2 --rnaseq_sets_dirs=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/rnaseq/ --BAMTOOLS_PATH=/home/s-amknut/GALBA/tools/bamtools/build/src/ --PROTHINT_PATH=/home/s-amknut/GALBA/tools/ProtHint/bin/ --prot_seq=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana/proteins/proteins.fa 




