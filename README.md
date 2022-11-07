Nothing is working, yet! Dummy code!

# GALBA User Guide

<u>Contact for Github Repository of GALBA at
https://github.com/Gaius-Augustus/GALBA:</u>

Katharina J. Hoff, University of Greifswald, Germany, katharina.hoff@uni-greifswald.de, +49 3834 420 4624


Author of GALBA
================

Katharina J. Hoff<sup name="aff1">[a, ](#aff1)</sup><sup name="aff2">[b](#aff2)</sup>

<b id="aff1">[a]</b> University of Greifswald, Institute for Mathematics and Computer Science, Walther-Rathenau-Str. 47, 17489 Greifswald, Germany

<b id="aff2">[b]</b> University of Greifswald, Center for Functional Genomics of Microbes, Felix-Hausdorff-Str. 8, 17489 Greifswald, Germany


Related Software
================

  * GALBA code was derived from BRAKER, a fully automated pipeline for predicting genes in the genomes of novel species with RNA-Seq data and a large-scale database of protein sequences (that must not necessarily be closely related to the target species) with GeneMark-ES/ET/EP/ETP and AUGUSTUS. BRAKER is available at https://github.com/Gaius-Augustus/BRAKER
  * TSEBRA can be used to combine GALBA gene sets with e.g. BRAKER gene sets. TSEBRA is available at https://github.com/Gaius-Augustus/TSEBRA .

Contents
========

-   [Authors](#authors-of-galba)
-   [What is GALBA?](#what-is-galba)
-   [Keys to successful gene prediction](#keys-to-successful-gene-prediction)
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
        -   [GALBA with proteins of short evolutionary distance](#galba-with-proteins-of-short-evolutionary-distance)
    -   [Description of selected GALBA command line options](#description-of-selected-galba-command-line-options)
        -   [--ab_initio](#--ab_initio)
        -   [--augustus_args=--some\_arg=bla](#--augustus_args--some_argbla)
        -   [--cores=INT](#--coresint)
        -   [--softmasking](#--softmasking)
        -   [--useexisting](#--useexisting)
        -   [--crf](#--crf)
        -   [--lambda=int](#--lambdaint)
	    -   [--makehub --email=your@mail.de](#--makehub---emailyourmailde)
-   [Output of GALBA](#output-of-galba)
-   [Example data](#example-data)
    -   [Data description](#data-description)
    -   [Testing GALBA with proteins of close homology](#testing-galba-with-proteins-of-close-homology)
    -   [Testing GALBA with pre-trained parameters](#testing-galba-with-pre-trained-parameters)
-   [Starting GALBA on the basis of previously existing GALBA runs](#starting-galba-on-the-basis-of-previously-existing-galba-runs)
-   [Bug reporting](#bug-reporting)
    -   [Reporting bugs on github](#reporting-bugs-on-github)
    -   [Common problems](#common-problems)
-   [Citing GALBA and software called by GALBA](#citing-galba-and-software-called-by-galba)
-   [License](#license)


What is GALBA?
===============

The rapidly growing number of sequenced genomes requires fully automated methods for accurate gene structure annotation. Here, we provide a fully automated gene pipeline that trains AUGUSTUS<sup name="a3">[R3, ](#f3)</sup><sup name="a4">[R4](#f4)</sup> for a novel species and subsequently predicts genes with AUGUSTUS in the genome of that species. GALBA uses protein sequences of a closely related species to generate a training gene set for AUGUSTUS with either miniprot<sup name="a1">[R1, ](#f1)</sup> or GenomeThreader<sup name="a2">[R2](#f2)</sup>. After training, GALBA uses the evidence from training genes during gene prediction.

<div class="panel panel-warning">
**Important**
{: .panel-heading}
<div class="panel-body">
Please note that the popular BRAKER<sup name="a5">[R5](#f5)</sup><sup name="a6">[R6](#f6)</sup> pipeline might produce more accurate results. Instead of using protein sequences of only one closely related species, BRAKER is capable of using proteins for a large sequence database where the species in the database must not necessarily be closely related to the target species. BRAKER can also incorporate RNA-Seq data. In contrast to GALBA, BRAKER achieves high gene prediction accuracy even in the absence of the annotation of very closely related species (and in the absence of RNA-Seq data).
</div>
</div>



Keys to successful gene prediction
==================================

-   Use a high quality genome assembly. If you have a huge number of very short scaffolds in your genome assembly, those short scaffolds will likely increase runtime dramatically but will not increase prediction accuracy.

-   Use simple scaffold names in the genome file (e.g. ```>contig1``` will work better than ```>contig1my custom species namesome putative function /more/information/  and lots of special characters %&!*(){}```). Make the scaffold names in all your fasta files simple before running any alignment program.

-   In order to predict genes accurately in a novel genome, the genome should be masked for repeats. This will avoid the prediction of false positive gene structures in repetitive and low complexitiy regions. In case of AUGUSTUS, softmasking (i.e. putting repeat regions into lower case letters and all other regions into upper case letters) leads to better results than hardmasking (i.e. replacing letters in repetitive regions by the letter `N` for unknown nucleotide). If the genome is softmasked, use the `--softmasking` flag of `galba.pl`.

-   Always check gene prediction results before further usage! You can e.g. use a genome browser for visual inspection of gene models in context with extrinsic evidence data. GALBA supports the generation of track data hubs for the UCSC Genome Browser with MakeHub for this purpose.

Overview running GALBA
====================================

GALBA mainly features semi-unsupervised, protein sequence evidence data supported training of AUGUSTUS with integration of extrinsic evidence in the final gene prediction step. GALBA can be used either with Miniprot or GenomeThreader as protein spliced aligner.


![galba2-sidetrack-b\[fig5\]](docs/figs/galba2_gth.png)

Figure a: training AUGUSTUS on the basis of spliced alignment information from proteins of a very closely related species against the target genome with Miniprot or GenomeThreader.


Installation
============

Supported software versions
---------------------------

At the time of release, this GALBA version was tested with:

-   AUGUSTUS 3.4.0 <sup name="g3">[F3, ](#g3)</sup><sup name="g4">[F4, ](#g4)</sup>

-   GenomeThreader 1.7.0<sup name="a2">[R2](#f2)</sup>

-   Miniprot<sup name="a1">[R1](#f1)</sup>

-   DIAMOND 0.9.24<sup name="a7">[R7](#f7)</sup>

-   cdbfasta 0.99

-   cdbyank 0.981

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

If you do not have root permissions on the Linux machine, try setting up an **Anaconda** (<https://www.anaconda.com/distribution/>) environment as follows:

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

Subsequently install GALBA and other software "as usual" while being in your conda environment.

### GALBA components

GALBA is a collection of Perl and Python scripts and a Perl module. The main script that will be called in order to run GALBA is `galba.pl`. Additional Perl and Python components are:

-   `align2hints.pl`

-   `filterIntronsFindStrand.pl`

-   `startAlign.pl`

-   `helpMod.pm`

-   `findGenesInIntrons.pl`

-   `downsample_traingenes.pl`

-   `ensure_n_training_genes.py`

All scripts (files ending with `*.pl` and `*.py`) that are part of GALBA must be executable in order to run GALBA. This should already be the case if you download GALBA from GitHub. Executability may be overwritten if you e.g. transfer GALBA on a USB-stick to another computer. In order to check whether required files are executable, run the following command in the directory that contains GALBA Perl scripts:

    ls -l *.pl *.py

The output should be similar to this:

```
    -rwxr-xr-x 1 katharina katharina  18191 Mai  7 10:25 align2hints.pl
    -rwxr-xr-x 1 katharina katharina   6090 Feb 19 09:35 galba_cleanup.pl
    -rwxr-xr-x 1 katharina katharina 408782 Aug 17 18:24 galba.pl
    -rwxr-xr-x 1 katharina katharina   5024 Mai  7 10:25 downsample_traingenes.pl
    -rwxr-xr-x 1 katharina katharina   5024 Mai  7 10:23 ensure_n_training_genes.py
    -rwxr-xr-x 1 katharina katharina   4542 Apr  3  2019 filter_augustus_gff.pl
    -rwxr-xr-x 1 katharina katharina   5754 Mai  7 10:25 filterIntronsFindStrand.pl
    -rwxr-xr-x 1 katharina katharina   7765 Mai  7 10:25 findGenesInIntrons.pl
    -rwxr-xr-x 1 katharina katharina   4679 Jan  9 13:55 merge_transcript_sets.pl
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

GALBA calls upon various bioinformatics software tools that are not part of GALBA. Some tools are obligatory, i.e. GALBA will not run at all if these tools are not present on your system. Other tools are optional. Please install all tools that are required for running GALBA in the mode of your choice.

### Mandatory tools

#### AUGUSTUS

Download AUGUSTUS from its master branch at <https://github.com/Gaius-Augustus/Augustus>. Unpack AUGUSTUS and install AUGUSTUS according to AUGUSTUS `README.TXT`. ***Do not use outdated AUGUSTUS versions from other sources, e.g. Debian package or Bioconda package! GALBA highly depends in particular on an up-to-date Augustus/scripts directory, and other sources are often lagging behind.***

You should compile AUGUSTUS on your own system in order to avoid problems with versions of libraries used by AUGUSTUS. Compilation instructions are provided in the AUGUSTUS `README.TXT` file (`Augustus/README.txt`).

AUGUSTUS consists of `augustus`, the gene prediction tool, additional C++ tools located in `Augustus/auxprogs` and Perl scripts located in `Augustus/scripts`. Perl scripts must be executable (see instructions in section [GALBA components](#executability). GALBA does not use any of the `auxprogs`. 

Since GALBA is a pipeline that trains AUGUSTUS, i.e. writes species specific parameter files, GALBA needs writing access to the configuration directory of AUGUSTUS that contains such files (`Augustus/config/`). If you install AUGUSTUS globally on your system, the `config` folder will typically not be writable by all users. Either make the directory where `config` resides recursively writable to users of AUGUSTUS, or copy the `config/` folder (recursively) to a location where users have writing permission.

AUGUSTUS will locate the `config` folder by looking for an environment variable `$AUGUSTUS_CONFIG_PATH`. If the `$AUGUSTUS_CONFIG_PATH`
environment variable is not set, then GALBA will look in the path `../config` relative to the directory in which it finds an AUGUSTUS
executable. Alternatively, you can supply the variable as a command line argument to GALBA (`--AUGUSTUS_CONFIG_PATH=/your_path_to_AUGUSTUS/Augustus/config/`). We recommend that you export the variable e.g. for your current bash
session:

```
    export AUGUSTUS_CONFIG_PATH=/your_path_to_AUGUSTUS/Augustus/config/
```

In order to make the variable available to all Bash sessions, add the above line to a startup script, e.g. `~/.bashrc`.

##### Important:

GALBA expects the entire `config` directory of AUGUSTUS at `$AUGUSTUS_CONFIG_PATH`, i.e. the subfolders `species` with its contents (at least `generic`) and `extrinsic`! Providing a writable but empty folder at `$AUGUSTUS_CONFIG_PATH` will not work for GALBA. If you need
to separate augustus binary and `$AUGUSTUS_CONFIG_PATH`, we recommend that you recursively copy the un-writable config contents to a writable location.

If you have a system-wide installation of AUGUSTUS at `/usr/bin/augustus`, an unwritable copy of `config` sits at `/usr/bin/augustus_config/`. The folder `/home/yours/` is writable to you. Copy with the following command (and additionally set the then required variables):

```
cp -r /usr/bin/Augustus/config/ /home/yours/
export AUGUSTUS_CONFIG_PATH=/home/yours/augustus_config
export AUGUSTUS_BIN_PATH=/usr/bin
export AUGUSTUS_SCRIPTS_PATH=/usr/bin/augustus_scripts
```

##### Modification of $PATH

Adding directories of AUGUSTUS binaries and scripts to your `$PATH` variable enables your system to locate these tools,
automatically. It is not a requirement for running GALBA to do this, because GALBA will try to guess them from the location of another
environment variable (`$AUGUSTUS_CONFIG_PATH`), or both directories can be supplied as command line arguments to `galba.pl`, but we recommend to add them to your `$PATH` variable. For your current bash session, type:

```
    PATH=:/your_path_to_augustus/bin/:/your_path_to_augustus/scripts/:$PATH
    export PATH
```

For all your BASH sessions, add the above lines to a startup script (e.g.`~/.bashrc`).

#### GenomeThreader

This tool is required, only, if you would like to run protein to genome alignments with GALBA using GenomeThreader. Download GenomeThreader from <http://genomethreader.org/>. Unpack and install according to `gth/README`.

GALBA will try to locate the GenomeThreader executable by using an environment variable `$ALIGNMENT_TOOL_PATH`. Alternatively, this can be supplied as command line argument (`--ALIGNMENT_TOOL_PATH=/your/path/to/gth`).

#### Miniprot


This tool is required, only, if you would like to run protein to genome alignments with GALBA using Miniprot. ... 


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


#### Tools from UCSC

If you switch `--UTR=on`, bamToWig.py will require the following tools that can be downloaded from <http://hgdownload.soe.ucsc.edu/admin/exe>:

 * twoBitInfo

 * faToTwoBit

It is optional to install these tools into your $PATH. If you don't, and you switch `--UTR=on`, bamToWig.py will automatically download them into the working directory.

#### MakeHub

If you wish to automaticaly generate a track data hub of your GALBA run, the MakeHub software, available at <https://github.com/Gaius-Augustus/MakeHub> is required. Download the software (either by running `git clone https://github.com/Gaius-Augustus/MakeHub.git`, or by
picking a release from <https://github.com/Gaius-Augustus/MakeHub/releases>. Extract the release package if you downloaded a release (e.g. `unzip MakeHub.zip` or `tar -zxvf MakeHub.tar.gz`.

GALBA will try to locate the make_hub.py script by using an environment variable `$MAKEHUB_PATH`. Alternatively, this can be supplied as command line argument (`--MAKEHUB_PATH=/your/path/to/MakeHub/`). GALBA can also try to guess the location of MakeHub on your system.

Running GALBA
===============

Different GALBA pipeline modes
--------------------------------

In the following, we describe the GALBA calls for Miniprot and GenomeThreader. In general, we recommend that you run GALBA on genomic sequences that have been softmasked for Repeats. If your genome has been softmasked, include the `--softmasking` flag in your GALBA call!

### GALBA with Miniprot

...

### GALBA with GenomeThreader

For running GALBA with GenomeThreader, type:

```
    galba.pl --species=yourSpecies --genome=genome.fasta \
       --prot_seq=proteins.fa --prg=gth \
       --ALIGNMENT_TOOL_PATH=/path/to/gth/binary \
       --trainFromGth
```

It is possible to generate protein alignments externally, prior running GALBA, itself. The compatible command for running GenomeThreader prior running GALBA, is:

```
    gth -genomic genome.fa  -protein protein.fa -gff3out \
       -skipalignmentout -o gth.aln
```

In order to use such externally created alignment files, run:

```
    galba.pl --species=yourSpecies --genome=genome.fasta \
       --prot_aln=proteins.aln --prg=gth --trainFromGth
```

It is also possible to run GALBA in this mode using an already prepared hints file. In this case, run:

```
    galba.pl --species=yourSpecies --genome=genome.fasta \
       --hints=hints.gff --prg=gth --trainFromGth
```

Format of the hints file should look like this:

```
    chrName   gth2h   CDSpart 105984  106633  .     -    .    src=P;grp=FBpp0285205;pri=4
    chrName   gth2h   start   106646  106648  .     -    .    src=P;grp=FBpp0285205;pri=4
```

Supported features in column 3 are intron, CDSpart, start, stop.


Description of selected GALBA command line options
----------------------------------------------------

Please run `galba.pl --help` to obtain a full list of options.

### --ab\_initio

Compute AUGUSTUS *ab initio* predictions in addition to AUGUSTUS predictions with hints (additional output files: `augustus.ab_initio.*`. This may be useful for estimating the quality of training gene parameters when inspecting predictions in a Browser.

### --augustus\_args="--some\_arg=bla"

One or several command line arguments to be passed to AUGUSTUS, if several arguments are given, separate them by whitespace, i.e. `"--first_arg=sth --second_arg=sth"`. This may be be useful if you know that gene prediction in your particular species benefits from a particular AUGUSTUS argument during the prediction step.

### --cores=INT

Specifies the maximum number of cores that can be used during computation. GALBA has to run some steps on a single core, others can take advantage of multiple cores. If you use more than 8 cores, this will not speed up all parallelized steps, in particular, the time consuming `optimize_augustus.pl` will not use more than 8 cores. However, if you don’t mind some cores being idle, using more than 8 cores will speed up other steps.

### --softmasking

Softmasking option for soft masked genome files. (Disabled by default.)

### --useexisting

Use the present config and parameter files if they exist for 'species'; will overwrite original parameters if GALBA performs an AUGUSTUS training.

### --crf

Execute CRF training for AUGUSTUS; resulting parameters are only kept for final predictions if they show higher accuracy than HMM parameters. This increases runtime!

### --lambda=int

Change the parameter $\lambda$ of the Poisson distribution that is used for downsampling training genes according to their number of introns (only genes with up to 5 introns are downsampled). The default value is $\lambda=2$. You might want to set it to 0 for organisms that mainly have single-exon genes. (Generally, single-exon genes contribute less value to increasing AUGUSTUS parameters compared to genes with many exons.)

### --makehub --email=your@mail.de

If `--makehub` and `--email=your@mail.de` (with your valid e-mail adress) are provided, a track data hub for visualizing results with the UCSC Genome Browser will be generated using MakeHub (<https://github.com/Gaius-Augustus/MakeHub>).


Output of GALBA
=================

GALBA produces several important output files in the working directory.

-   augustus.hints.gtf: Genes predicted by AUGUSTUS with hints from given extrinsic evidence. This file will be missing if GALBA was run with the option `--esmode`.

-   augustus.ab_initio.gtf: Genes predicted by AUGUSTUS in *ab initio* mode in GTF-format. The file will always be present if AUGUSTUS has been run with the option `--esmode`. Otherwise, it will only be present if GALBA was run with the option `--AUGUSTUS_ab_initio`.

-   hintsfile.gff: The extrinsic evidence data extracted from RNAseq.bam and/or protein data.

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

If the `--makehub` option was used and MakeHub is available on your system, a hub directory beginning with the name `hub_` will be created. Copy this directory to a publicly accessible web server. A file `hub.txt` resides in the directory. Provide the link to that file to the UCSC Genome Browser for visualizing results.

Example data
============

An example data set is contained in the directory `GALBA/example`.

In case you have trouble accessing that file, there's also a copy available from another server:

The example data set was not compiled in order to achieve optimal prediction accuracy, but in order to quickly test pipeline components. The small subset of the genome used in these test examples is not long enough for GALBA training to work well.

Data description
----------------

Data corresponds to the last 1,000,000 nucleotides of _Arabidopsis thaliana_'s chromosome Chr5, split into 8 artificial contigs.

The protein sequences are a subset of [OrthoDB v10 plants proteins](https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz).

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

Testing GALBA with GenomeThreader
---------------------------------

    galba.pl --genome genome.fa --prot_seq proteins.fa --prg gth \
        --trainFromGth --softmasking --cores N


This test is implemented in `test1.sh`, expected runtime is ~7 minutes. The fast runtime of this test is mostly caused by generating a low number of training genes. Note that this approach does not scale well with increasing genome size and the number of proteins in a protein database.


Testing GALBA with pre-trained parameters
------------------------------------------

The training step of all pipelines can be skipped with the option `--skipAllTraining`. This means, only AUGUSTUS predictions will be performed, using pre-trained, already existing parameters. For example, you can predict genes with the command:

```
    galba.pl --genome=genome.fa --bam RNAseq.bam --species=arabidopsis \
        --skipAllTraining --softmasking --cores N
```

This test is implemented in `test3.sh`, expected runtime is ~1 minute.

Bug reporting
=============

Before reporting bugs, please check that you are using the most recent versions of AUGUSTUS and GALBA. Also, check the list of [Common problems](#common-problems), and the Issue list on GitHub before reporting bugs. We do monitor open issues on GitHub. Sometimes, we are unable to help you, immediately, but we try hard to solve your problems.

Reporting bugs on GitHub
------------------------

If you found a bug, please open an issue at <https://github.com/Gaius-Augustus/GALBA/issues> (or contact katharina.hoff@uni-greifswald.de).

Information worth mentioning in your bug report:

Check in `galba/yourSpecies/galba.log` at which step `galba.pl` crashed.

There are a number of other files that might be of interest, depending on where in the pipeline the problem occurred. Some of the following files will not be present if they did not contain any errors.

-   `galba/yourSpecies/hintsfile.gff` - is this file empty? If yes, something went wrong during hints generation - does this file contain hints from source “b2h” and of type “intron”?

-   `galba/yourSpecies/align_gth/*err` - errors reported by the alignment tools gth

-   `galba/yourSpecies/genbank.good.gb` - try a “grep -c LOCUS genbank.good.gb” to determine the number of training genes for training AUGUSTUS, should not be low

-   `galba/yourSpecies/errors/firstetraining.stderr` - contains errors from first iteration of training AUGUSTUS

-   `galba/yourSpecies/errors/secondetraining.stderr` - contains errors from second iteration of training AUGUSTUS

-   `galba/yourSpecies/errors/optimize_augustus.stderr` - contains errors optimize\_augustus.pl (additional training set for AUGUSTUS)

-   `galba/yourSpecies/errors/augustus*.stderr` - contain AUGUSTUS execution errors

-   `galba/yourSpecies/startAlign.stderr` - if you provided a protein fasta file and `--prg` option and this file is not empty, something went wrong during protein alignment

-   `galba/yourSpecies/startAlign.stdout` - may give clues on at which point protein alignment went wrong

Common problems
---------------

-   *There are duplicate Loci in the `train.gb` file (after using GenomeThreader)!*

    This issue arises if outdated versions of AUGUSTUS and GALBA are used. Solution: Please update AUGUSTUS and GALBA from github (<https://github.com/Gaius-Augustus/Augustus>, <https://github.com/Gaius-Augustus/GALBA>).

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

    -   Stanke. M., Schöffmann, O., Morgenstern, B. and Waack, S. (2006). Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources. BMC Bioinformatics 7, 62.

    - Buchfink, B., Xie, C., Huson, D.H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods 12:59-60.

-   If GALBA was executed with Miniprot, cite:

    -   ...

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

<b id="f1">[R1]</b> Li, H. (2022) “Protein-to-genome alignment with miniprot.” arXiv:2210.08052v1.[↩](#a1)

<b id="f2">[R2]</b> Gremme, G. 2013. “Computational Gene Structure Prediction.” PhD thesis, Universität Hamburg.[↩](#a2)

<b id="f3">[R3]</b> Stanke, M., Schöffmann, O., Morgenstern, B., and Waack., S. 2006. “Gene Prediction in Eukaryotes with a Generalized Hidden Markov Model That Uses Hints from External Sources.” *BMC Bioinformatics* 7 (1). BioMed Central: 62.[↩](#a3)

<b id="f4"><R4></b> Hoff, K. and Stanke, M. 2019. “Predicting genes in single genomes with AUGUSTUS.“ *Current Protocols in Bioinformatics*, 65(1), e57.[↩](#a4)

<b id="f5">[R5]</b> Bruna, T., Hoff, K. J.,  Lomsadze, A., Stanke, M., and Borodvsky, M. 2021. “BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database." *NAR Genomics and Bioinformatics* 3(1):lqaa108.[↩](#a5)

<b id="f6">[R6]</b> Hoff, K. J., Lange, S., Lomsadze, A., Borodovsky, M., and Stanke, M. 2015. “BRAKER1: Unsupervised Rna-Seq-Based Genome Annotation with Genemark-et and Augustus.” *Bioinformatics* 32 (5), 767--69.[↩](#a6)

<b id="f7">[R7]</b>Buchfink, B., Xie, C., and Huson, D. H. 2015. Fast and sensitive protein alignment using DIAMOND. *Nature Methods*, 12(1), 59-60.[↩](#a7)

<b id="f8">[R8]</b> Hoff, K.J. 2019. MakeHub: Fully automated generation of UCSC Genome Browser Assembly Hubs. *Genomics, Proteomics and Bioinformatics*, in press, preprint on bioarXive, doi: <https://doi.org/10.1101/550145>.[↩](#a8)
