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
	"getAnnoFasta.augustus.hints.stdout", "secondtest.stdout",
	"train.gb", 
	"aug_hints.lst", "aa2nonred.stdout", "augustus.hints.tmp.gtf", "fourthtest.stdout", "gbFilterEtraining.stdout",
	"genes.gtf", "genes_in_gb.gtf",
	"singlecds.hints", "stops.and.starts.gff", "train.gb.test", "train.gb.train", "train.gb.train.test",
	"train.gb.train.train", "traingenes.good.fa", "downsample_traingenes.log", "firstetraining.stdout", 
	"secondetraining.stdout", "startAlign_gth.log", "protein_alignment_gth.gff3", "ex1.cfg", "getAnnoFastaFromJoingenes.augustus.hints.stdout", 
	"genome.fa.cidx", "getAnnoFastaFromJoingenes.augustus.hints_tmp.stdout", "getAnnoFastaFromJoingenes.augustus.ab_initio_tmp.stdout",
	"augustus.ab_initio.tmp.gtf", "augustus.ab_initio.gff", "augustus.hints.tmp.gtf",
	"getAnnoFastaFromJoingenes.augustus.hints_hints.stdout", "getAnnoFastaFromJoingenes.augustus.ab_initio_.stdout",
	"getAnnoFastaFromJoingenes.augustus.hints_.stdout", "startAlign.stdout", "cmd.log", "etrain.bad.lst", "gene_stat.yaml",
	"good_genes.lst", "nonred.loci.lst", "nuc.fasta", "proteins.fa", 
    "train.f.gb", "traingenes.good.gtf", "traingenes.good.nr.fa", "uniqueSeeds.gtf", "genome.mpi", "pygustus_hints.out", "pygustus_hints.py",
	"gff2gbSmallDNA.stderr", "miniprot_representatives.gff", "protein_alignment_miniprot.aln", "protein_alignment_miniprot.gff", "hc.gff", 
	"miniprothint.gff", "thirdtest.stdout");


foreach(@files){
	if(-e $wdir."/".$_){
		print "Deleting file ".$wdir."/".$_."\n";
		unlink($wdir."/".$_);
	}
}

my $gth_index = 0;
while(-d $wdir."/align_gth".$gth_index){
	print "Deleting directory ".$wdir."/align_gth".$gth_index."\n";
	rmtree($wdir."/align_gth".$gth_index);
	$gth_index = $gth_index + 1;
}

# create a new directory for archived files
my $archivedir = $wdir."/archive";
if(not(-d $archivedir)){
	mkdir($archivedir);
}

# move files that are not essential output to archive directory, 
# there files are: augustus.hints.gff, train2.gb.train.test, train2.gb.train.train, train2.gb.test, train2.gb.train, train2.gb, miniprot_representatives.gtf, miniprot.gff
my @archivefiles = ("augustus.hints.gff", "train2.gb.train.test", "train2.gb.train.train", "train2.gb.test", "train2.gb.train", "train2.gb", "miniprot_representatives.gtf", "miniprot.gff");
foreach(@archivefiles){
	if(-e $wdir."/".$_){
		print "Moving file ".$wdir."/".$_." to archive directory\n";
		rename($wdir."/".$_, $archivedir."/".$_);
	}
}