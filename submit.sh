#!/bin/bash
#SBATCH -J my_job                            # job name
#SBATCH -N 1                                  # number of nodes
#SBATCH -n 16                                  # number of MPI processes, here 1 per node
#SBATCH --partition=snowball                 # choose nodes from partition
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

#Run pregalba.py for A. thaliana in mixed mode
WD=/home/nas-hs/projs/galba-isoseq/data/Arabidopsis_thaliana
singularity exec -B $WD galba.sif /path/to/pregalba.py -y /path/to/arabidopsis.yaml -t 16 --mixed 

#Run pregalba.py for M. musculus in isoseq mode
WD=/home/nas-hs/projs/galba-isoseq/data/Mus_musculus
singularity exec -B $WD galba.sif /path/to/pregalba.py -y /path/to/musmusculus.yaml -t 16 --isoseq 

