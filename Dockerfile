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

#galba
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/GALBA.git  && \
    git clone --depth=1 https://github.com/Gaius-Augustus/BRAKER.git && \
    cp  BRAKER/scripts/compute_accuracies.sh GALBA/scripts/compute_accuracies.sh && \
    cp BRAKER/scripts/compare_intervals_exact.pl GALBA/scripts/compare_intervals_exact.pl

FROM $BASE_CONTAINER

USER root

COPY --from=base /opt/ /opt/

ENV PATH=${PATH}:/opt/seqstats:/opt/cdbfasta:/opt/hisat2:/opt/diamond:/opt/TSEBRA/bin:/opt/MakeHub:/opt/miniprot:/opt/GALBA/scripts:/opt/miniprot-boundary-scorer:/opt/miniprothint

ENV AUGUSTUS_CONFIG_PATH=/usr/share/augustus/config/

# augustus
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

# patch Augustus scripts (because Debian package is often outdated, this way we never need to worry)
RUN cd /tmp/ && mkdir removeme && cd removeme && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/optimize_augustus.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/aa2nonred.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/gff2gbSmallDNA.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/new_species.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/filterGenesIn.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/filterGenes.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/join_mult_hints.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/randomSplit.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/getAnnoFastaFromJoingenes.py && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/gtf2gff.pl && \
    chmod a+x *.pl && \ 
    cd /usr/share/augustus/scripts && \
    cp /tmp/removeme/* . && \
    rm -rf /tmp/removeme

ENV AUGUSTUS_BIN_PATH=/usr/bin/
ENV AUGUSTUS_SCRIPTS_PATH=/usr/share/augustus/scripts/

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

WORKDIR "${HOME}"
