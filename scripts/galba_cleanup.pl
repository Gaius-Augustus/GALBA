#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# galba_cleanup.pl                                                                                 #
# Script that deletes all files produced by GALBA that are not deleted by the galba.pl script,     #
# itself, and that are usually not required for downstream analysis, unless you wish to debug      #
# something.                                                                                       #
#                                                                                                  #
# Author: Katharina Hoff                                                                           #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
####################################################################################################

use Getopt::Long;
use File::Path;
use File::Spec::Functions qw(rel2abs);
use strict;
use warnings;

my $usage = <<'ENDUSAGE';

galba_cleanup.pl      delete galba.pl output files that are usually not required for 
                      downstream analysis

SYNOPSIS

galba_cleanup.pl --wdir=WDIR

OPTIONS

--wdir=WDIR         output directory of galba.pl job

--help	            Display this help message

ENDUSAGE

my $help;
my $wdir;
my $homedir = $ENV{HOME};

GetOptions(
	'wdir=s'  => \$wdir,
    'help!' => \$help
);

if ($help) {
	print $usage;
    exit(0);
}

if(not(defined($wdir))){
	print "ERROR: in file " . __FILE__ ." at line "
            . __LINE__ . "\n" . "No Working directory provided!"
            . "(option --wdir=WDIR)\n";
    print $usage;
	exit(1)
}

if($wdir =~ m/^~/){
	$wdir =~ s/^~//;
	$wdir = File::Spec->catfile($homedir, $wdir);
}

if(not($wdir =~ m/^\//)){
	$wdir = rel2abs($wdir);
}

if(not(-d $wdir)){
	my @t = split(/\//, $wdir);
	$wdir = "";
	for(my $i=0; $i<(scalar(@t) - 1); $i++){
		$wdir .= "/".$t[$i]
	}
	if(not(-d $wdir)){
		print "ERROR: in file " . __FILE__ ." at line "
            . __LINE__ . "\n" . "wdir $wdir is not a directory!\n";
	exit(1)
	}
}

my @files = ("firsttest.stdout", "genome.fa", "getAnnoFasta.augustus.ab_initio.stdout", 
	"getAnnoFasta.augustus.hints.stdout", "optimize_augustus.stdout", "secondtest.stdout",
	"train.gb", 
	"aug_hints.lst", "aa2nonred.stdout", "augustus.hints.tmp.gtf", "bothutr.lst", "fourthtest.stdout", "gbFilterEtraining.stdout",
	"genes.gtf", "genes_in_gb.gtf", "gff2gbSmallDNA.utr.stdout", "hints.job.lst", "merged.bam", "merged.s.bam", "merged.wig", 
	"randomSplit_utr1.log", "randomSplit_utr2.log", "rnaseq2utr.stdout", "rnaseq.utr.hints", "samtools_sort_before_wig.stdout",
	"singlecds.hints", "stops.and.starts.gff", "train.gb.test", "train.gb.train", "train.gb.train.test",
	"train.gb.train.train", "traingenes.good.fa", "traingenes.gtf", "utr.aa2nonred.stdout",
	"utr.gb", "utr.gb.test", "utr.gb.train", "utr.gb.train.test", "utr.gb.train.train", "utr_genes_in_gb.fa",
	"utr_genes_in_gb.nr.fa", "utr.nonred.loci.lst", "utrs.gff", "downsample_traingenes.log", "firstetraining.stdout", 
	"secondetraining.stdout", "startAlign_gth.log", "protein_alignment_gth.gff3", "ex1.cfg", "getAnnoFastaFromJoingenes.augustus.hints.stdout", 
	"genome.fa.cidx", "getAnnoFastaFromJoingenes.augustus.hints_tmp.stdout", "getAnnoFastaFromJoingenes.augustus.ab_initio_tmp.stdout",
	"augustus.ab_initio.tmp.gtf", "augustus.ab_initio.gff", "augustus.hints.gff", "augustus.hints.tmp.gtf",
	"getAnnoFastaFromJoingenes.augustus.hints_hints.stdout", "getAnnoFastaFromJoingenes.augustus.ab_initio_.stdout",
	"getAnnoFastaFromJoingenes.augustus.hints_.stdout", "startAlign.stdout", "augustus.hints_iter1.aa", "augustus.hints_iter1.codingseq",
	"augustus.hints_iter1.gff", "augustus.tmp1.gff", "augustus.tmp2.gff", "cmd.log", "etrain.bad.lst", "gene_stat.yaml",
	"good_genes.lst", "hintsfile_iter1.gff", "nonred.loci.lst", "nuc.fasta", "prevHints.gff", "proteins.fa", "seed_proteins.faa",
        "train.f.gb", "traingenes.good.gtf", "traingenes.good.nr.fa", "uniqueSeeds.gtf", "genome.mpi", "protein_alignment_miniprot.gff");

my @dirs = ("genome_split", "Spaln", "augustus_tmp", "augustus_files_hints", "diamond");

foreach(@files){
	if(-e $wdir."/".$_){
		print "Deleting file ".$wdir."/".$_."\n";
		unlink($wdir."/".$_);
	}
}

foreach(@dirs){
	if(-d $wdir."/".$_){
		print "Deleting directory ".$wdir."/".$_."\n";
		rmtree($wdir."/".$_);
	}
}

my $gth_index = 0;
while(-d $wdir."/align_gth".$gth_index){
	print "Deleting directory ".$wdir."/align_gth".$gth_index."\n";
	rmtree($wdir."/align_gth".$gth_index);
	$gth_index = $gth_index + 1;
}

