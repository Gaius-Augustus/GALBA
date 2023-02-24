#!/bin/bash

# Author: Katharina J. hoff
# Contact: katharina.hoff@uni-greifswald.de
# Date: Jan 12th 2023

# Copy this script into the folder where you want to execute it, e.g.:
# singularity exec -B $PWD:$PWD galba.sif cp /opt/GALBA/example/singularity-tests/test1.sh .
# Then run "bash test1.sh".

# Check whether galba.sif is available
if [[ -z "${GALBA_SIF}" ]]; then
    echo ""
    echo "Variable GALBA_SIF is undefined."
    echo "First, build the sif-file with \"singularity build galba.sif docker://katharinahoff/galba-notebook:latest\""
    echo ""
    echo "After building, export the GALBA_SIF environment variable on the host as follows:"
    echo ""
    echo "export GALBA_SIF=\$PWD/galba.sif"
    echo ""
    echo "You will have to modify the export statement if galba.sif does not reside in \$PWD."
    echo ""
    exit 1
fi

# Check whether singularity exists
if ! command -v singularity &> /dev/null
then
    echo "Singularity could not be found."
    echo "On some HPC systems you can load it with \"module load singularity\"."
    echo "If that fails, please install singularity."
    echo "Possibly you misunderstood how to run this script. Before running it, please copy it to the directory where you want to execute it by e.g.:"
    echo "singularity exec -B \$PWD:\$PWD galba.sif cp /opt/GALBA/example/singularity-tests/test1.sh ."
    echo "Then execute on the host with \"bash test1.sh\"".
    exit 1
fi

# remove output directory if it already exists
wd=test1
if [ -d $wd ]; then
    rm -r $wd
fi

singularity exec -B ${PWD}:${PWD} ${GALBA_SIF} galba.pl --genome=/opt/GALBA/example/genome.fa --prot_seq=/opt/GALBA/example/proteins.fa --workingdir=${wd} --threads 8 --skipOptimize
            # Important: the option --skipOptimize should never be applied to a real life run!!!
            # It as only introduced to speed up the test. Please delete it from the script if you use it for real data analysis.
