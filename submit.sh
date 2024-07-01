#!/bin/bash
#SBATCH -J my_job                            # job name
#SBATCH -N 1                                  # number of nodes
#SBATCH -n 1                                   # number of MPI processes, here 1 per node
#SBATCH --partition=snowball         # choose nodes from partition
#SBATCH -o %j.out                            # stdout file name (%j: job ID)
#SBATCH -e %j.err                             # stderr file name (%j: job ID)
#SBATCH -t 01:00:00                        # max run time (hh:mm:ss), max 72h!
#SBATCH --mail-type=end
#SBATCH --mail-user=s-amknut@uni-greifswald.de

## optional environment variables
echo "On which nodes it executes:"
echo $SLURM_JOB_NODELIST
echo " "
echo "jobname: $SLURM_JOB_NAME"

## load modules
#module load gcc/5.2.0   # for instance

./bin/main.py /home/s-amknut/GALBA/bin/testGenome.fa /home/s-amknut/GALBA/bin/testGenome /home/s-amknut/GALBA/bin/testReads.fa