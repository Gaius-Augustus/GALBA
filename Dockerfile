# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
ARG OWNER=jupyter
ARG BASE_CONTAINER=$OWNER/minimal-notebook
FROM $BASE_CONTAINER as base

# Fix: https://github.com/hadolint/hadolint/wiki/DL4006
# Fix: https://github.com/koalaman/shellcheck/wiki/SC3014
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root

RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN apt update && \
    apt-get install -y --no-install-recommends \
    man-db \
    g++ \
    less \
    zlib1g-dev \
    && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# out of my way
RUN rm -rf /opt && mkdir /opt

RUN cd /opt && \ 
    git clone --recursive https://github.com/clwgg/seqstats && \
    cd seqstats && \
    make 

# cdbfasta
RUN cd /opt && \ 
    git clone --depth=1 https://github.com/gpertea/cdbfasta.git && \
    cd cdbfasta && \
    make 

# diamond
RUN cd /opt && \
    mkdir diamond && \
    cd diamond && \
    wget http://github.com/bbuchfink/diamond/releases/download/v2.0.15/diamond-linux64.tar.gz && \
    tar -xf diamond-linux64.tar.gz && \
    rm diamond-linux64.tar.gz

# tsebra
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/TSEBRA && \
    cd TSEBRA && \
    git checkout v1.0.3

# makehub
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/MakeHub.git && \
    cd MakeHub && \
    git checkout braker3 && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedToBigBed && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/genePredCheck && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/faToTwoBit && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/gtfToGenePred && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/hgGcPercent && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/ixIxx && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/twoBitInfo && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/genePredToBed && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/genePredToBigGenePred && \
    chmod u+x bedToBigBed genePredCheck faToTwoBit gtfToGenePred hgGcPercent ixIxx  twoBitInfo wigToBigWig genePredToBed genePredToBigGenePred make_hub.py

#miniprot
RUN cd /opt && \
    git clone --depth=1 https://github.com/lh3/miniprot && \
    cd miniprot && \
    make

#miniprot-boundary-scorer
RUN cd /opt && \
    git clone --depth=1 https://github.com/tomasbruna/miniprot-boundary-scorer.git && \
    cd miniprot-boundary-scorer && \
    make

#miniprothint
RUN cd /opt && \
    git clone --depth=1 https://github.com/tomasbruna/miniprothint.git

# get AUGUSTUS compilation dependencies
# Install required packages
RUN apt-get update
RUN apt-get install -y build-essential wget git autoconf

# Install dependencies for AUGUSTUS comparative gene prediction mode (CGP)
RUN apt-get install -y libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev
RUN apt-get install -y libsqlite3-dev libmysql++-dev

# Install dependencies for the optional support of gzip compressed input files
RUN apt-get install -y libboost-iostreams-dev zlib1g-dev

# Install dependencies for bam2hints and filterBam
RUN apt-get install -y libbamtools-dev

# Install additional dependencies for bam2wig
RUN apt-get install -y samtools libhts-dev

# Install additional dependencies for homGeneMapping and utrrnaseq
RUN apt-get install -y libboost-all-dev

# compile augustus from source because of segmentation fault
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/Augustus.git && \
    cd Augustus && \
    make && \
    cd scripts && \
    chmod a+x *.pl && \
    chmod a+x *.py


FROM $BASE_CONTAINER

USER root

COPY --from=base /opt/ /opt/

ENV PATH=${PATH}:/opt/seqstats:/opt/cdbfasta:/opt/hisat2:/opt/diamond:/opt/TSEBRA/bin:/opt/MakeHub:/opt/miniprot:/opt/GALBA/scripts:/opt/miniprot-boundary-scorer:/opt/miniprothint

# AUGUSTUS does need several libraries that are now gone, re-install them:
RUN apt-get update --yes && \
    apt-get install -y libboost-iostreams-dev zlib1g-dev libboost-all-dev libboost-all-dev libbamtools-dev

ENV AUGUSTUS_CONFIG_PATH=/opt/Augustus/config/

# augustus, install only in order to get the dependencies, will uninstall augustus later on
RUN apt update && \ 
    apt install -yq  augustus augustus-data augustus-doc \
    # for latex labels
    cm-super \
    dvipng \
    # for matplotlib anim
    ffmpeg \
    time && \
    apt clean all && \
    fix-permissions "${AUGUSTUS_CONFIG_PATH}"

# perl & dependencies
RUn apt update && \
    apt install -yq libyaml-perl \
                    libhash-merge-perl \
                    libparallel-forkmanager-perl \
                    libscalar-util-numeric-perl \
                    libclass-data-inheritable-perl \
                    libexception-class-perl \
                    libtest-pod-perl \
                    libfile-which-perl \
                    libmce-perl \
                    libthread-queue-perl \
                    libmath-utils-perl \
                    libscalar-list-utils-perl && \
    apt clean all

USER ${NB_UID}

# only python installations can be done as a normal user
RUN mamba install --quiet -c bioconda -c anaconda --yes \
    biopython && \
    mamba clean  --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

# Pygustus in Pip is at version 0.8.3
RUN pip install pygustus && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

USER root

RUN apt-get remove -y augustus augustus-data augustus-doc

ENV AUGUSTUS_BIN_PATH=/opt/Augustus/bin/
ENV AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts/
ENV PATH=${PATH}:/opt/Augustus/scripts/:/opt/Augustus/bin/

#galba
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/GALBA.git  && \
    git clone --depth=1 https://github.com/Gaius-Augustus/BRAKER.git && \
    cp BRAKER/scripts/compute_accuracies.sh GALBA/scripts/compute_accuracies.sh && \
    cp BRAKER/scripts/compare_intervals_exact.pl GALBA/scripts/compare_intervals_exact.pl

USER ${NB_UID}

WORKDIR "${HOME}"
