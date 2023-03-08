# GALBA User Guide

![Docker Pulls](https://img.shields.io/docker/pulls/katharinahoff/galba-notebook)

:warning: If you pulled from dockerhub between 24th and 28th of February, you may want to upgrade your container to the latest version, as soon as possible. There were a couple of serious bugs with the Pygustus integration... apologies. It has been fixed.

<u>Contact for Github Repository of GALBA at
https://github.com/Gaius-Augustus/GALBA:</u>

Katharina J. Hoff, University of Greifswald, Germany, katharina.hoff@uni-greifswald.de, +49 3834 420 4624


Authors of GALBA
================

Katharina J. Hoff<sup name="aff1">[a, ](#aff1)</sup><sup name="aff2">[b](#aff2)</sup>, Heng Li<sup name="aff3">[c, ](#aff3)</sup><sup name="aff4">[d ](#aff4)</sup>, Tomas Bruna<sup name="aff5">[e](#aff1)</sup>, Lars Gabriel<sup name="aff1">[a, ](#aff1)</sup><sup name="aff2">[b](#aff2)</sup>, and Mario Stanke<sup name="aff1">[a, ](#aff1)</sup><sup name="aff2">[b](#aff2)</sup>

<b id="aff1">[a]</b> University of Greifswald, Institute for Mathematics and Computer Science, Walther-Rathenau-Str. 47, 17489 Greifswald, Germany

<b id="aff2">[b]</b> University of Greifswald, Center for Functional Genomics of Microbes, Felix-Hausdorff-Str. 8, 17489 Greifswald, Germany

<b id="aff3">[c]</b> Dana-Farber Cancer Institute, 450 Brookline Ave, Boston, MA 02215, USA

<b id="aff4">[d]</b> Harvard Medical School, 10 Shattuck St, Boston, MA 02215, USA

<b id="aff5">[e]</b> Joint Genome Institute, Lawrence Berkeley National Laboratory, USA

Acknowledgements
===============

GALBA code was derived from the BRAKER code, where a similar pipeline for using GenomeThreader with BRAKER was once published in <sup name="a9">[R9](#f9)</sup>. We hereby acknowledge the contributions of all BRAKER authors to the code that GALBA was derived from, and we are grateful for funding for BRAKER development by the National Institutes of Health (NIH) grant GM128145, which indirectly also supported development of GALBA.


Related Software
================

  * GALBA code was derived from BRAKER, a fully automated pipeline for predicting genes in the genomes of novel species with RNA-Seq data and a large-scale database of protein sequences (that must not necessarily be closely related to the target species) with GeneMark-ES/ET/EP/ETP and AUGUSTUS. BRAKER is available at https://github.com/Gaius-Augustus/BRAKER
  * TSEBRA can be used to combine GALBA gene sets with BRAKER gene sets. TSEBRA is available at https://github.com/Gaius-Augustus/TSEBRA .

Contents
========

-   [Authors](#authors-of-galba)
-   [Acknowledgements](#acknowledgements)
-   [What is GALBA?](#what-is-galba)
-   [Keys to successful gene prediction](#keys-to-successful-gene-prediction)
-   [Singularity Image](#singularity-image)
-   [Installation](#installation)
    -   [Supported software versions](#supported-software-versions)
    -   [GALBA](#galba)
        -   [Perl pipeline dependencies](#perl-pipeline-dependencies)
        -   [GALBA components](#galba-components)
        -   [Bioinformatics software dependencies](#bioinformatics-software-dependencies)
            -   [Mandatory tools](#mandatory-tools)
            -   [Optional tools](#optional-tools)
-   [Running GALBA](#running-galba)
    -   [GALBA pipeline modes](#different-galba-pipeline-modes)
    -   [Description of selected GALBA command line options](#description-of-selected-galba-command-line-options)
        -   [--AUGUSTUS_ab_initio](#--AUGUSTUS_ab_initio)
        -   [--augustus_args=--some\_arg=bla](#--augustus_args--some_argbla)
        -   [--threads=INT](#--threadsint)
        -   [--crf](#--crf)
        -   [--lambda=int](#--lambdaint)
	    -   [--makehub --email=your@mail.de](#--makehub---emailyourmailde)
-   [Output of GALBA](#output-of-galba)
-   [Example data](#example-data)
    -   [Data description](#data-description)
    -   [Testing GALBA with Miniprot](#testing-galba-with-miniprot)
    -   [Testing GALBA with GenomeThreader](#testing-galba-with-genomethreader)
    -   [Testing GALBA with pre-trained parameters](#testing-galba-with-pre-trained-parameters)
-   [Accuracy](#accuracy)
-   [Bug reporting](#bug-reporting)
    -   [Reporting bugs on github](#reporting-bugs-on-github)
    -   [Common problems](#common-problems)
-   [Citing GALBA and software called by GALBA](#citing-galba-and-software-called-by-galba)
-   [License](#license)


What is GALBA?
===============

The rapidly growing number of sequenced genomes requires fully automated methods for accurate gene structure annotation. Here, we provide a fully automated gene pipeline that trains AUGUSTUS<sup name="a3">[R3, ](#f3)</sup><sup name="a4">[R4](#f4)</sup> for a novel species and subsequently predicts genes with AUGUSTUS in the genome of that species. GALBA uses the protein sequences of **one closely related species** to generate a training gene set for AUGUSTUS with either miniprot<sup name="a1">[R1, ](#f1)</sup> or GenomeThreader<sup name="a2">[R2](#f2)</sup>. After training, GALBA uses the evidence from protein to genome alignment during gene prediction. 

:warning: Please note that the popular BRAKER<sup name="a5">[R5, ](#f5)</sup><sup name="a6">[R6](#f6)</sup> pipeline will very likely produce more accurate results than GALBA. Instead of using protein sequences of only one closely related species, **BRAKER is capable of using proteins from a large sequence database** where the species in the database must not necessarily be closely related to the target species. BRAKER can also incorporate RNA-Seq data. In contrast to GALBA, BRAKER achieves high gene prediction accuracy even in the absence of the annotation of very closely related species (and in the absence of RNA-Seq data). Before deciding to use GALBA, please read the [Accuracy](#accuracy) section.

**If you are not sure which pipeline to use: GALBA or BRAKER? The answer is: use BRAKER, first!**

GALBA is named after Servius Sulpicius Galba, who ruled the Roman Empire only for a short time, before he was murdered. The named seems appropriate, because both BRAKER2 and also the soon published BRAKER3 achieve higher accuracy than GALBA ever will.

Keys to successful gene prediction
==================================

-   Use a high quality genome assembly. If you have a huge number of very short scaffolds in your genome assembly, those short scaffolds will likely increase runtime dramatically but will not increase prediction accuracy.

-   Use simple scaffold names in the genome file (e.g. ```>contig1``` will work better than ```>contig1my custom species namesome putative function /more/information/  and lots of special characters %&!*(){}```). Make the scaffold names in all your fasta files simple before running any alignment program.

-   In order to predict genes accurately in a novel genome, the genome should be masked for repeats. This will avoid the prediction of false positive gene structures in repetitive and low complexitiy regions. In the case of AUGUSTUS, softmasking (i.e., putting repeat regions into lower case letters and all other regions into upper case letters) leads to better results than hardmasking (i.e., replacing letters in repetitive regions by the letter `N` for unknown nucleotide). GALBA always treats genomes as softmasked for repeats!

-   Always check gene prediction results before further usage! You can, e.g. use a genome browser for visual inspection of gene models in context with extrinsic evidence data. GALBA supports the generation of track data hubs for the UCSC Genome Browser with MakeHub for this purpose.

Overview running GALBA
====================================

GALBA mainly features semi-unsupervised, protein sequence evidence data supported training of AUGUSTUS with integration of extrinsic evidence in the final gene prediction step. GALBA can be used either with Miniprot or GenomeThreader as a protein spliced aligner. Miniprot is our preferred aligner because it continues to undergo development and is faster than GenomeThreader.

![galba-miniprot\[fig1\]](docs/figs/galba_miniprot.png)

Figure a: training AUGUSTUS on the basis of spliced alignment information from proteins of a very closely related species against the target genome with miniprot.

![galba-gth\[fig2\]](docs/figs/galba_gth.png)

Figure b:training AUGUSTUS on the basis of spliced alignment information from proteins of a very closely related species against the target genome with GenomeThreader.

Singularity Image
=================

The easiest way to run GALBA is using singuarlity. We provide a docker container to build a singularity image (tested with singularity version 3.10.0-dirty). **We only include Miniprot in Docker & Singularity!** GenomeThreader is not included.

Build as follows (requires 1.2 GB disk space):

```
singularity build galba.sif docker://katharinahoff/galba-notebook:latest
```

Execute GALBA from galba.sif like this (i.e. it automatically mounts the user's home directory on the host system):

```
singularity exec galba.sif galba.pl
```

Running GALBA in Singularity outside of $HOME
---------------------------------------------

If you want to execute galba.sif while mounting a different directory, e.g. mounting $PWD, then you need to be aware of the following: GALBA needs a writable `$AUGUSTUS_CONFIG_PATH` environment variable. By default, the `AUGUSTUS_CONFIG_PATH` is `/usr/share/augustus/config` in the sif container, which is not writable. Therefore, GALBA attempts to automatically copy the contents of `/usr/share/augustus/config` into $HOME/.augustus on the host system. If this fails, it copies to $PWD.

If you want to re-use AUGUSTUS parameters trained by GALBA in a later run with --skipAllTraining, you must either mount the same $PWD as during training, or you must manually copy the $PWD/.augustus to the location that you will mount for the second run.

What's included in docker/singuarlity GALBA
--------------------------------------------

Among others, the containers include the following software that is useful in context of working with GALBA:

   * galba.pl
   * augustus
   * TSEBRA
   * make_hub.py & UCSC dependencies
   * miniprot
   * diamond
   * cdbfasta
   * seqstats

Installation
============

The long way around is to manually install all dependencies of GALBA.

Supported software versions
---------------------------

At the time of release, this GALBA version was tested with:

-   AUGUSTUS 3.4.0 <sup name="g3">[R3, ](#g3)</sup><sup name="g4">[R4](#g4)</sup>

-   Pygustus v0.8.0-alpha

-   Miniprot 0.5-r181<sup name="a1">[R1](#f1)</sup>

-   GenomeThreader 1.7.3<sup name="a2">[R2](#f2)</sup>

-   DIAMOND 0.9.230<sup name="a7">[R7](#f7)</sup>

-   cdbfasta 0.99

-   cdbyank 0.981

-   miniprot-boundary-scorer

-   Miniprothint

GALBA
-------

### Perl pipeline dependencies

Running GALBA requires a Linux-system with `bash` and Perl. Furthermore, GALBA requires the following CPAN-Perl modules to be
installed:

-   `File::Spec::Functions`

-   `Hash::Merge`

-   `List::Util`

-   `MCE::Mutex`

-   `Module::Load::Conditional`

-   `Parallel::ForkManager`

-   `POSIX`

-   `Scalar::Util::Numeric`

-   `YAML`

-   `Math::Utils`

-   `File::HomeDir`


On Ubuntu, for example, install the modules with CPANminus<sup name="g4">[F4](#g4)</sup>: `sudo cpanm Module::Name`, e.g. `sudo cpanm Hash::Merge`.

GALBA also uses a Perl module `helpMod.pm` that is not available on CPAN. This module is part of the GALBA release and does not require separate installation.

If you do not have root permissions on the Linux machine, try running it from the [Singularity image](#singularity-image), or try setting up an **Anaconda** (<https://www.anaconda.com/distribution/>) environment as follows:

```
wget https://repo.anaconda.com/archive/Anaconda3-2018.12-Linux-x86_64.sh
bash bin/Anaconda3-2018.12-Linux-x86_64.sh # do not install VS (needs root privileges)
conda install -c anaconda perl
conda install -c bioconda perl-app-cpanminus
conda install -c bioconda perl-hash-merge
conda install -c bioconda perl-parallel-forkmanager
conda install -c bioconda perl-scalar-util-numeric
conda install -c bioconda perl-yaml
conda install -c bioconda perl-class-data-inheritable
conda install -c bioconda perl-exception-class
conda install -c bioconda perl-test-pod
conda install -c anaconda biopython
conda install -c bioconda perl-file-which # skip if you are not comparing to reference annotation
conda install -c bioconda perl-mce
conda install -c bioconda perl-list-util
conda install -c bioconda perl-math-utils
conda install -c bioconda cdbtools
```

Subsequently, install GALBA and other software "as usual" while in your conda environment.

### GALBA components

GALBA is a collection of Perl and Python scripts and a Perl module. The main script that will be called in order to run GALBA is `galba.pl`. Additional Perl and Python components are:

-   `aln2hints.pl`

-   `filterIntronsFindStrand.pl`

-   `startAlign.pl`

-   `helpMod.pm`

-   `downsample_traingenes.pl`

-   `galba_cleanup.pl`

All scripts (files ending with `*.pl` and `*.py`) that are part of GALBA must be executable in order to run GALBA. This should already be the case if you download GALBA from GitHub. Executability may be overwritten if you, e.g. transfer GALBA on a USB-stick to another computer. In order to check whether required files are executable, run the following command in the directory that contains GALBA Perl scripts:

    ls -l *.pl *.py

The output should be similar to this:

```
    -rwxr-xr-x 1 katharina katharina  18191 Mai  7 10:25 aln2hints.pl
    -rwxr-xr-x 1 katharina katharina   6090 Feb 19 09:35 galba_cleanup.pl
    -rwxr-xr-x 1 katharina katharina 408782 Aug 17 18:24 galba.pl
    -rwxr-xr-x 1 katharina katharina   5024 Mai  7 10:25 downsample_traingenes.pl
    -rwxr-xr-x 1 katharina katharina   5754 Mai  7 10:25 filterIntronsFindStrand.pl
    -rwxr-xr-x 1 katharina katharina  41674 Mai  7 10:25 startAlign.pl
```

It is important that the `x` in `-rwxr-xr-x` is present for each script. If that is not the case, run

    `chmod a+x *.pl *.py`

in order to change file attributes.

You may find it helpful to add the directory in which GALBA perl scripts reside to your `$PATH` environment variable. For a single bash session, enter:

```
    PATH=/your_path_to_galba/:$PATH
    export PATH
```

To make this `$PATH` modification available to all bash sessions, add the above lines to a startup script (e.g.`~/.bashrc`).

Bioinformatics software dependencies
------------------------------------

GALBA calls upon various bioinformatics software tools that are not part of GALBA. Some tools are mandatory, i.e. GALBA will not run at all if these tools are not present on your system. Other tools are optional. Please install all tools that are required for running GALBA in the mode of your choice.

### Mandatory tools

#### AUGUSTUS

Download AUGUSTUS from its master branch at <https://github.com/Gaius-Augustus/Augustus>. Unpack AUGUSTUS and install AUGUSTUS according to AUGUSTUS `README.TXT`. ***Do not use outdated AUGUSTUS versions from other sources, e.g. Debian package or the Bioconda package! GALBA highly depends in particular on an up-to-date Augustus/scripts directory, and other sources are often lagging behind.***

You should compile AUGUSTUS on your own system in order to avoid problems with versions of libraries used by AUGUSTUS. Compilation instructions are provided in the AUGUSTUS `README.TXT` file (`Augustus/README.txt`).

AUGUSTUS consists of `augustus`, the gene prediction tool, additional C++ tools located in `Augustus/auxprogs` and Perl scripts located in `Augustus/scripts`. Perl scripts must be executable (see instructions in section [GALBA components](#executability). GALBA does not use any of the `auxprogs`.

Since GALBA is a pipeline that trains AUGUSTUS, i.e., writes species specific parameter files, GALBA needs write access to the configuration directory of AUGUSTUS that contains such files (`Augustus/config/`). If you install AUGUSTUS globally on your system, the `config` folder will typically not be writable by all users. Either make the directory where `config` resides recursively writable to users of AUGUSTUS, or copy the `config/` folder (recursively) to a location where users have writing permission.

AUGUSTUS will locate the `config` folder by looking for the environment variable `$AUGUSTUS_CONFIG_PATH`. If the `$AUGUSTUS_CONFIG_PATH`
environment variable is not set, then GALBA will look in the path `../config` relative to the directory in which it finds an AUGUSTUS
executable. Alternatively, you can supply the variable as a command line argument to GALBA (`--AUGUSTUS_CONFIG_PATH=/your_path_to_AUGUSTUS/Augustus/config/`). We recommend that you export the variable, e.g., for your current bash
session:

```
    export AUGUSTUS_CONFIG_PATH=/your_path_to_AUGUSTUS/Augustus/config/
```

In order to make the variable available to all Bash sessions, add the above line to a startup script, e.g., `~/.bashrc`.

##### Important:

GALBA expects the entire `config` directory of AUGUSTUS at `$AUGUSTUS_CONFIG_PATH`, i.e. the subfolders `species` with its contents (at least `generic`) and `extrinsic`! Providing a writable but empty folder at `$AUGUSTUS_CONFIG_PATH` will not work for GALBA. If you need
to separate augustus binary and `$AUGUSTUS_CONFIG_PATH`, we recommend that you recursively copy the un-writable config contents to a writable location.

If you have a system-wide installation of AUGUSTUS at `/usr/bin/augustus`, an unwritable copy of `config` sits at `/usr/bin/augustus_config/`. The folder `${HOME}` is writable to you. Copy with the following command (and additionally set the then required variables):

```
cp -r /usr/bin/Augustus/config/ ${HOME}/augustus_config
export AUGUSTUS_CONFIG_PATH=${HOME}/augustus_config
export AUGUSTUS_BIN_PATH=/usr/bin
export AUGUSTUS_SCRIPTS_PATH=/usr/bin/augustus_scripts
```
GALBA automatically sets the environment variables if run from the Singuarity image.

##### Modification of $PATH

Adding directories of AUGUSTUS binaries and scripts to your `$PATH` variable enables your system to locate these tools,
automatically. It is not a requirement for running GALBA to do this, because GALBA will try to guess them from the location of another
environment variable (`$AUGUSTUS_CONFIG_PATH`), or both directories can be supplied as command line arguments to `galba.pl`, but we recommend to add them to your `$PATH` variable. For your current bash session, type:

```
PATH=:/your_path_to_augustus/bin/:/your_path_to_augustus/scripts/:$PATH
export PATH
```

For all your BASH sessions, add the above lines to a startup script (e.g.`~/.bashrc`).

#### Pygustus

This tool is required for parallelization of Augustus. Install Pygustus (<https://github.com/Gaius-Augustus/pygustus>) with pip (or pip3):

```
pip install pygustus
```

Important: this version of GALBA relies on pygustus 0.8.3-alpha. It will not be compatible with older versions.

#### Miniprot

This tool is only required, if you would like to run protein to genome alignments with GALBA using Miniprot. Download Miniprot from <https://github.com/lh3/miniprot>:

```
git clone https://github.com/lh3/miniprot.git
cd miniprot
make
```

GALBA will try to locate the Miniprot executable by using an environment variable `$MINIPROT_PATH`. Alternatively, this can be supplied as command line argument (`--MINIPROT_PATH=/your/path/to/miniprot/`).

#### miniprot-boundary-scorer

This tool is only required, if you would like to run protein to genome alignments with GALBA using Miniprot. Download miniprot-boundary-scorer from <https://github.com/tomasbruna/miniprot-boundary-scorer>:

```
git clone https://github.com/tomasbruna/miniprot-boundary-scorer.git
cd miniprot-boundary-scorer
make
```

GALBA will try to locate the miniprot-boundary-scorer executable by using an environment variable `$SCORER_PATH`. Alternatively, this can be supplied as command line argument (`--SCORER_PATH=/your/path/to/miniprot-boudary-scorer/`).

#### miniprotprothint

This tool is only required, if you would like to run protein to genome alignments with GALBA using Miniprot. Download miniprothint from <https://github.com/tomasbruna/miniprothint>:

```
git clone https://github.com/tomasbruna/miniprothint.git
```

GALBA will try to locate the miniprothint.py executable by using an environment variable `$MINIPROTHINT_PATH`. Alternatively, this can be supplied as command line argument (`--MINIPROTHINT_PATH=/your/path/to/miniprothint/`).

#### GenomeThreader

This tool is only required, if you would like to run protein to genome alignments with GALBA using GenomeThreader. Download GenomeThreader from <http://genomethreader.org/>. Unpack and install according to `gth/README`.

GALBA will try to locate the GenomeThreader executable by using an environment variable `$GENOMETHREADER_PATH`. Alternatively, this can be supplied as command line argument (`--GENOMETHREADER_PATH=/your/path/to/gth/`).

Please be aware that miniprot achieves higher accuracy than GenomeThreader in GALBA! See section on [Accuracy](#accuracy).


#### Python3

On Ubuntu, Python3 is usually installed by default, `python3` will be in your `$PATH` variable, by default, and GALBA will automatically locate it. However, you have the option to specify the `python3` binary location in two other ways:

1.  Export an environment variable `$PYTHON3_PATH`, e.g. in your `~/.bashrc` file:

        export PYTHON3_PATH=/path/to/python3/

2.  Specify the command line option `--PYTHON3_PATH=/path/to/python3/` to `galba.pl`.


#### DIAMOND

DIAMOND is used for removal of redundant training genes.

Obtain and unpack DIAMOND as follows:

```
    wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz
    tar xzf diamond-linux64.tar.gz
```

If already in your `$PATH` variable, GALBA will find diamond, automatically. Otherwise, GALBA can locate the diamond binary either by using an environment variable `$DIAMOND_PATH`, or by taking a command line argument (`--DIAMOND_PATH=/your_path_to_diamond`). In order to set the environment variable e.g. for your current bash session, type:

```
    export DIAMOND_PATH=/your_path_to_diamond/
```

Add the above line to a startup script (e.g. `~/.bashrc`) in order to set the environment variable for all bash sessions.

### Optional tools

#### Biopython

If Biopython is installed, GALBA can generate FASTA-files with coding sequences and protein sequences predicted by AUGUSTUS and generate track data hubs for visualization of a GALBA run with MakeHub <sup name="a8">[R8](#f8)</sup>.
These are optional steps. The first can be disabled with the command-line flag `--skipGetAnnoFromFasta`, the second can be activated by using the command-line options `--makehub --email=your@mail.de`, Biopython is not required if neither of these optional steps shall be performed.

On Ubuntu, install Python3 package manager with:

    `sudo apt-get install python3-pip`

Then, install Biopython with:

    `sudo pip3 install biopython`

#### cdbfasta

cdbfasta and cdbyank are required by GALBA for correcting AUGUSTUS genes with in frame stop codons (spliced stop codons) using the AUGUSTUS script fix_in_frame_stop_codon_genes.py. This can be skipped with `--skip_fixing_broken_genes`.

On Ubuntu, install cdbfasta with:

    `sudo apt-get install cdbfasta`

For other systems, you can for example obtain cdbfasta from <https://github.com/gpertea/cdbfasta>, e.g.:

```
        git clone https://github.com/gpertea/cdbfasta.git`
        cd cdbfasta
        make all
```

On Ubuntu, cdbfasta and cdbyank will be in your `$PATH` variable after installation, and GALBA will automatically locate them. However, you have the option to specify the `cdbfasta` and `cdbyank` binary location in two other ways:

1.  Export an environment variable `$CDBTOOLS_PATH`, e.g. in your `~/.bashrc` file:

```
        export CDBTOOLS_PATH=/path/to/cdbtools/
```

2.  Specify the command line option `--CDBTOOLS_PATH=/path/to/cdbtools/` to `galba.pl`.

#### MakeHub

If you wish to automaticaly generate a track data hub of your GALBA run, the MakeHub software, available at <https://github.com/Gaius-Augustus/MakeHub> is required. Download the software (either by running `git clone https://github.com/Gaius-Augustus/MakeHub.git`, or by
picking a release from <https://github.com/Gaius-Augustus/MakeHub/releases>. Extract the release package if you downloaded a release (e.g. `unzip MakeHub.zip` or `tar -zxvf MakeHub.tar.gz`.

GALBA will try to locate the make_hub.py script by using an environment variable `$MAKEHUB_PATH`. Alternatively, this can be supplied as command line argument (`--MAKEHUB_PATH=/your/path/to/MakeHub/`). GALBA can also try to guess the location of MakeHub on your system.

System dependencies
-------------------
galba.pl uses getconf to see how many threads can be run on your system. On Ubuntu, you can install it with:

```
sudo apt-get install libc-bin
```

Running GALBA
===============

In the following, we describe the GALBA calls for Miniprot and GenomeThreader. In general, we recommend that you run GALBA on genomic sequences that have been softmasked for Repeats.

### GALBA with Miniprot

For running GALBA with Miniprot, type:

```
    galba.pl --species=yourSpecies --genome=genome.fasta \
       --prot_seq=proteins.fa
```

### GALBA with GenomeThreader

For running GALBA with GenomeThreader, type:

```
    galba.pl --species=yourSpecies --genome=genome.fasta \
       --prot_seq=proteins.fa --prg=gth
```

Description of selected GALBA command line options
----------------------------------------------------

Please run `galba.pl --help` to obtain a full list of options.

### --prg

Use either `miniprot` (default) or `gth` (for GenomeThreader) to generate training genes and hints.

### --AUGUSTUS_ab_initio

Compute AUGUSTUS *ab initio* predictions in addition to AUGUSTUS predictions with hints (additional output files: `augustus.ab_initio.*`. This may be useful for estimating the quality of training gene parameters when inspecting predictions in a Browser.

### --augustus_args="--some_arg=bla"

One or several command line arguments to be passed to AUGUSTUS, if several arguments are given, separate them by whitespace, i.e. `"--first_arg=sth --second_arg=sth"`. This may be be useful if you know that gene prediction in your particular species benefits from a particular AUGUSTUS argument during the prediction step.

### --threads=INT

Specifies the maximum number of threads that can be used during computation. GALBA has to run some steps on a single thread, others can take advantage of multiple threads. If you use more than 8 threads, this will not speed up all parallelized steps, in particular, the time consuming `optimize_augustus.pl` will not use more than 8 threads.

### --crf

Execute CRF training for AUGUSTUS; resulting parameters are only kept for final predictions if they show higher accuracy than HMM parameters. This increases runtime!

### --lambda=int

Change the parameter $\lambda$ of the Poisson distribution that is used for downsampling training genes according to their number of introns (only genes with up to 5 introns are downsampled). The default value is $\lambda=2$. You might want to set it to 0 for organisms that mainly have single-exon genes. (Generally, single-exon genes contribute less value to increasing AUGUSTUS parameters compared to genes with many exons.)

### --makehub --email=your@mail.de

If `--makehub` and `--email=your@mail.de` (with your valid e-mail adress) are provided, a track data hub for visualizing results with the UCSC Genome Browser will be generated using MakeHub (<https://github.com/Gaius-Augustus/MakeHub>).


Output of GALBA
=================

GALBA produces several important output files in the working directory.

-   augustus.hints.gtf: Genes predicted by AUGUSTUS with hints from given extrinsic evidence.

-   augustus.ab_initio.gtf: Genes predicted by AUGUSTUS in *ab initio* mode in GTF-format. The file will always be present if AUGUSTUS has been run with the option `--esmode`. Otherwise, it will only be present if GALBA was run with the option `--AUGUSTUS_ab_initio`.

-   hintsfile.gff: The extrinsic evidence data extracted from protein data.

AUGUSTUS output files may be present with the following name endings and formats:

-   GTF-format is always produced.

-   GFF3-format is produced if the flag `--gff3` was specified to GALBA.

-   Coding sequences in FASTA-format are produced if the flag `--skipGetAnnoFromFasta` was not set.

-   Protein sequence files in FASTA-format are produced if the flag `--skipGetAnnoFromFasta` was not set.

For details about gtf format, see <http://www.sanger.ac.uk/Software/formats/GFF/>. A GTF-format file contains one line per predicted exon. Example:

```
    HS04636 AUGUSTUS initial   966 1017 . + 0 transcript_id "g1.1"; gene_id "g1";
    HS04636 AUGUSTUS internal 1818 1934 . + 2 transcript_id "g1.1"; gene_id "g1";
```

The columns (fields) contain:

```
    seqname source feature start end score strand frame transcript ID and gene ID
```

If the `--makehub` option was used and MakeHub is available on your system, a hub directory beginning with the name `hub_` will be created. Copy this directory to a publicly accessible web server. A file `hub.txt` resides in the directory. Provide the link to that file to the UCSC Genome Browser for visualizing results. MakeHub is included in the Singularity image.

Example data
============

An example data set is contained in the directory `GALBA/example`.

The example data set was not compiled in order to achieve optimal prediction accuracy, but in order to quickly test pipeline components. The small subset of the genome used in these test examples is not long enough for GALBA training to work well.

Data description
----------------

Data corresponds to the last 1,000,000 nucleotides of _Arabidopsis thaliana_'s chromosome Chr5, split into 8 artificial contigs.

The protein sequences are a subset from *Arabidopsis lyrata subsp. lyrata* genome assembly GCA_000004255.1.

List of files:

-   `genome.fa` - genome file in fasta format
-   `proteins.fa` - protein sequences in fasta format

The below given commands assume that you configured all paths to tools by exporting bash variables or that you have the necessary tools in your $PATH.

The example data set also contains scripts `tests/test*.sh` that will execute below listed commands for testing GALBA with the example data set. You find example results of AUGUSTUS in the folder `results/test*`. Be aware that GALBA contains several parts where random variables are used, i.e. results that you obtain when running the tests may not be exactly identical. To compare your test results with the reference ones, you can use the [compare_intervals_exact.pl](https://github.com/Gaius-Augustus/BRAKER/blob/master/scripts/compare_intervals_exact.pl) script from BRAKER as follows:

    # Compare CDS features
    compare_intervals_exact.pl --f1 augustus.hints.gtf --f2 ../../results/test${N}/augustus.hints.gtf --verbose
    # Compare transcripts
    compare_intervals_exact.pl --f1 augustus.hints.gtf --f2 ../../results/test${N}/augustus.hints.gtf --trans --verbose

We give runtime estimations derived from computing on *Intel(R) Xeon(R) CPU E5530 @ 2.40GHz*.

Testing GALBA with Miniprot
---------------------------

```
galba.pl --genome=genome.fa --prot_seq=proteins.fa \
    --skipOptimize --threads 8
```

This test is implemented in `test1.sh`, expected runtime is ~2 minutes. The fast runtime of this test is mostly caused by generating a low number of training genes, and by skipping an optimization step for AUGUSTUS training.

Testing GALBA with GenomeThreader
---------------------------------

```
galba.pl --genome=genome.fa --prot_seq=proteins.fa \
    --skipOptimize --prg=gth --workingdir=$wd \
    --threads 8
```

This test is implemented in `test2.sh`, expected runtime is ~2:15 minutes. The fast runtime of this test is mostly caused by generating a low number of training genes, and by skipping an optimization step for AUGUSTUS training.

Testing GALBA with pre-trained parameters
------------------------------------------

The training step of all pipelines can be skipped with the option `--skipAllTraining`. This means, only AUGUSTUS predictions will be performed, using pre-trained, already existing parameters. For example, you can predict genes with the command:

```
galba.pl --genome=genome.fa --prot_seq=proteins.fa \
    --skipAllTraining \
    --threads 8 --species=arabidopsis
```

This test is implemented in `test3.sh`, expected runtime is 2:30 minutes.

Accuracy
========

![galba-miniprot-fly\[fig3\]](docs/figs/galba_miniprot_fly.png)

Figure c: accuracy results of GALBA and BRAKER2 in *Drosophila melanogaster*. Shown are results from running GALBA with protein input from *D. simulans* (dsim), *D. erecta* (dere), *D. ananassae* (dana), *D. speudoobscura* (dpse), *D. willistoni* (dwil), *D. virilis* (dvir), and *D. grimshawi* (dgri). In addition, GALBA was run with proteins of the house fly (mdom). For "combo", the input proteins were a concatenation of dana, dpse, dwil, dere, and dgri. BRAKER2 results were computed with the OrthoDB Arthropoda partition, excluding *D. melanogaster* proteins.

:warning: As Figure c demonstrates, BRAKER2 (<https://github.com/Gaius-Augustus/BRAKER>) generates more accurate results than GALBA in fruit fly. We have computed similar results for *Arabidopsis thaliana* and *Caenorhabditis elegans*. Therefore, we generally recommend that you use BRAKER instead of GALBA!

There may be some special cases where GALBA obtains better results than BRAKER. For example, if you observe a "split gene" problem with BRAKER, and if you have proteins of a very close relative at hand, then GALBA may improve over BRAKER. Also, if you miss genes that are known in relatives of your species of interest in a BRAKER output, it might be worth trying GALBA and combining the resulting gene set and hintsfile with TSEBRA with a BRAKER output.

Accuracy of GALBA with GenomeThreader is generally lower than with miniprot. In many cases, GenomeThreader generates too few training genes for training AUGUSTUS, therefore we show only a subset of the data from Figure c in Figure d:

![galba-miniprot-gth-fly\[fig4\]](docs/figs/galba_miniprot_gth_fly.png)

Bug reporting
=============

Before reporting bugs, please check that you are using the most recent versions of AUGUSTUS and GALBA. Also, check the list of [Common problems](#common-problems), and the Issue list on GitHub before reporting bugs. We do monitor open issues on GitHub. Sometimes, we are unable to help you, immediately, but we try hard to solve your problems.

Reporting bugs on GitHub
------------------------

If you found a bug, please open an issue at <https://github.com/Gaius-Augustus/GALBA/issues> (or contact katharina.hoff@uni-greifswald.de).

Information worth mentioning in your bug report:

Check in `${wd}/galba.log` at which step `galba.pl` crashed.

There are a number of other files that might be of interest, depending on where in the pipeline the problem occurred. Some of the following files will not be present if they did not contain any errors.

-   `${wd}/hintsfile.gff` - is this file empty? If yes, something went wrong during hints generation - does this file contain hints from source “b2h” and of type “intron”?

-   `${wd}/genbank.good.gb` - try a “grep -c LOCUS genbank.good.gb” to determine the number of training genes for training AUGUSTUS, should not be low

-   `${wd}/errors/firstetraining.stderr` - contains errors from first iteration of training AUGUSTUS

-   `${wd}/errors/secondetraining.stderr` - contains errors from second iteration of training AUGUSTUS

-   `${wd}/errors/optimize_augustus.stderr` - contains errors optimize\_augustus.pl (additional training set for AUGUSTUS)

-   `${wd}/errors/augustus*.stderr` - contain AUGUSTUS execution errors

-   `${wd}/startAlign.stderr` - if you provided a protein fasta file and `--prg=gth` option and this file is not empty, something went wrong during protein alignment

-   `${wd}/startAlign.stdout` - may give clues on at which point protein alignment went wrong

Common problems
---------------

-   *[something] failed to execute!*

    When providing paths to software to GALBA, please use absolute, non-abbreviated paths. For example, GALBA might have problems with `--DIAMOND_PATH=./diamond/` or `--DIAMOND_PATH=~/diamond/`. Please use `DIAMOND_PATH=/full/absolute/path/to/diamond/`, instead. This applies to all path specifications as command line options to `galba.pl`. Relative paths and absolute paths will not pose problems if you export a bash variable, instead, or if you append the location of tools to your $PATH variable.

-   *GALBA cannot find the Augustus script XYZ...*

    Update Augustus from github with `git clone https://github.com/Gaius-Augustus/Augustus.git`. Do not use Augustus from other sources. GALBA is highly dependent on an up-to-date Augustus. Augustus releases happen rather rarely, updates to the Augustus scripts folder occur rather frequently.

-   *Does GALBA depend on Python3?*

    It does. The python scripts employed by GALBA are not compatible with Python2.

-   *Why does GALBA predict more genes than I expected?*

    If transposable elements (or similar) have not been masked appropriately, AUGUSTUS tends to predict those elements as protein coding genes. This can lead to a huge number genes. You can check whether this is the case for your project by BLASTing (or DIAMONDing) the predicted protein sequences against themselves (all vs. all) and counting how many of the proteins have a high number of high quality matches. You can use the output of this analysis to divide your gene set into two groups: the protein coding genes that you want to find and the repetitive elements that were additionally predicted.


Citing GALBA and software called by GALBA
=============================================

Since GALBA is a pipeline that calls several Bioinformatics tools, publication of results obtained by GALBA requires that not only GALBA is referred to, but also the tools that are called by GALBA. GALBA will output a file `what-to-cite.txt` in the GALBA working directory, informing you about which exact sources apply to your run.

-   Always cite:

    -   Hoff, K.J., Lomsadze, A., Borodovsky, M. and Stanke, M. (2019). Whole-Genome Annotation with BRAKER. Methods Mol Biol. 1962:65-95, doi: 10.1007/978-1-4939-9173-0_5.

    -   Hoff, K. and Stanke, M. 2019. “Predicting genes in single genomes with AUGUSTUS.“ *Current Protocols in Bioinformatics*, 65(1), e57.

    -   Stanke. M., Schöffmann, O., Morgenstern, B. and Waack, S. (2006). Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources. BMC Bioinformatics 7, 62.

    -   Buchfink, B., Xie, C., Huson, D.H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods 12:59-60.

-   If GALBA was executed with Miniprot, cite:

    -   Li, H. (2022) “Protein-to-genome alignment with miniprot.” arXiv:2210.08052v1.

-   If GALBA was executed with GenomeThreader, cite:

    -   Gremme, G. (2013). Computational Gene Structure Prediction. PhD thesis, Universität Hamburg.

-   If GALBA called MakeHub for creating a track data hub for visualization of GALBA results with the UCSC Genome Browser, cite:

    -   Hoff, K.J. (2019) MakeHub: Fully automated generation of UCSC Genome Browser Assembly Hubs. Genomics, Proteomics and Bioinformatics, in press 2020, preprint on bioarXive, doi: <https://doi.org/10.1101/550145>.


License
=======

All source code, i.e. `scripts/*.pl` or `scripts/*.py` are under the Artistic License (see <http://www.opensource.org/licenses/artistic-license.php>).

Footnotes
=========

<b id="g2">[F2]</b> Please use the latest version from the master branch of AUGUSTUS distributed by the original developers, it is available from github at <https://github.com/Gaius-Augustus/Augustus>. Problems have been reported from users that tried to run GALBA with AUGUSTUS releases maintained by third parties, i.e. Bioconda. [↩](#g2)

<b id="g4">[F4]</b> install with `sudo apt-get install cpanminus` [↩](#g4)


References
==========

<b id="f1">[R1]</b> Li, H. (2022) “Protein-to-genome alignment with miniprot.” arXiv:2210.08052v1.[↩](#a1)

<b id="f2">[R2]</b> Gremme, G. 2013. “Computational Gene Structure Prediction.” PhD thesis, Universität Hamburg.[↩](#a2)

<b id="f3">[R3]</b> Stanke, M., Schöffmann, O., Morgenstern, B., and Waack., S. 2006. “Gene Prediction in Eukaryotes with a Generalized Hidden Markov Model That Uses Hints from External Sources.” *BMC Bioinformatics* 7 (1). BioMed Central: 62.[↩](#a3)

<b id="f4">[R4]</b> Hoff, K. and Stanke, M. 2019. “Predicting genes in single genomes with AUGUSTUS.“ *Current Protocols in Bioinformatics*, 65(1), e57.[↩](#a4)

<b id="f5">[R5]</b> Bruna, T., Hoff, K. J.,  Lomsadze, A., Stanke, M., and Borodvsky, M. 2021. “BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database." *NAR Genomics and Bioinformatics* 3(1):lqaa108.[↩](#a5)

<b id="f6">[R6]</b> Hoff, K. J., Lange, S., Lomsadze, A., Borodovsky, M., and Stanke, M. 2015. “BRAKER1: Unsupervised Rna-Seq-Based Genome Annotation with Genemark-et and Augustus.” *Bioinformatics* 32 (5), 767--69.[↩](#a6)

<b id="f7">[R7]</b>Buchfink, B., Xie, C., and Huson, D. H. 2015. Fast and sensitive protein alignment using DIAMOND. *Nature Methods*, 12(1), 59-60.[↩](#a7)

<b id="f8">[R8]</b> Hoff, K.J. 2019. MakeHub: Fully automated generation of UCSC Genome Browser Assembly Hubs. *Genomics, Proteomics and Bioinformatics*, in press, preprint on bioarXive, doi: <https://doi.org/10.1101/550145>.[↩](#a8)

<b id="f9">[R9]</b>Hoff, K.J., Lomsadze, A., Borodovsky, M. and Stanke, M. (2019). Whole-Genome Annotation with BRAKER. Methods Mol Biol. 1962:65-95, doi: 10.1007/978-1-4939-9173-0_5.[↩](#a9)
