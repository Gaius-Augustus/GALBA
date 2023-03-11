#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# aln2hints.pl -   generate hints from spaln [O0 (=gff3)], exonerate, genomeThreader (gth)         #
#                  or scipio output                                                                #
#                                                                                                  #
# Script is derived from align2hints.pl in BRAKER.                                                 #
#                                                                                                  #
# Authors: Katharina Hoff                                                                          #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Last modification: November 4th 2022                                                             #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
# Usage:                                                                                           #
# align2hints.pl [OPTIONS] --in=align.gff3 --out=hintsfile.gff \                                   #
#                --prg=miniprot|gth                                                                #
#                                                                                                  #
####################################################################################################

use strict;
use warnings;
use Getopt::Long;
use File::Spec::Functions qw(rel2abs);
use Cwd;

my $usage = <<'ENDUSAGE';

aln2hints.pl    generate hints from miniprot or 
                  GenomeThreader (gth) output
                  Miniprot run like this:
                      miniprot -ut16 --gtf genome.fa proteins.fa > miniprot.gtf
                  GenomeThreader run like this: 
                      gth -genomic genome.fa  -protein protein.fa -gff3out \
                         -skipalignmentout ... -o gthfile

SYNOPSIS

align2hints.pl [OPTIONS] --in=align.gff --out=hintsfile.gff \
                         --prg=miniprot|gth

  --in                 input file from gth (gff3), miniprot (gtf)
  --out                contains CDSpart, CDS and intron hints


OPTIONS

    --help                   Print this help message.
    --CDSpart_cutoff=n       This many bp are cut off of each CDSpart hint
                             w.r.t. the cds (default 15).
    --maxintronlen=n         Alignments with longer gaps are discarded
                             (default 350000).
    --minintronlen=n         Alignments with gaps shorter than this and longer
                             than maxgaplen are discarded (default 41).
    --priority=n             Priority of hint group (default 4).
    --prg=s                  Alignment program of input file, either 'gth' or
                             'miniprot'.
    --source=s               Source identifier (default 'P')
    --CDS                    Do not output CDSpart hints, but complete CDS
                             hints.
    --genome_file=s          if prg is miniprot and start hints shall be
                             created, the genome file from which the
                             alignments were generated, must be specified.
    --version                print version of aln2hints.pl

Format:
  seqname <TAB> source <TAB> feature <TAB> start <TAB> end <TAB> score <TAB>
     strand <TAB> frame <TAB> src=source;grp=target_protein;pri=priority


DESCRIPTION

  Example:

    aln2hints.pl [OPTIONS] --in=align.gtf --out=hintsfile.gff --prg=miniprot

ENDUSAGE

my $version = 1.0;    # version of aln2hints.pl
my $printVersion;
my $alignfile;        # alignment input file
my $CDSpart_cutoff = 15;           # cutoff for CDSpart hints
my $CDSpartid      = "CDSpart";    # abbreviate to decrease file size
my $dir;               # working superdirectory where programme is called from
my $hintsfilename;     # hints file output name
my $intron_end;        # end of intron (only for gth and spaln)
my $intron_score;      # intron score
my $intron_start;      # start of intron (only for gth and spaln)
my $intron_threshold;  # Threshold for current programme
my $intron_threshold_gth
    = 0.7;    # introns from gth with a score below this are discarded
my $maxintronlen = 350000;    # maximal intron length
my $minintronlen = 41;        # default minimal intron length
my $parent;                   # current parent
my $prevParent = "noP"; # previous parent
my $prevScore  = 0;     # previous exon/CDS score for calculating intron score
# positions in query protein (scipio only), to determine intron for scipio
my ( $qstart, $qend, $prevQend );
my $prgsrc;    # source programme (exonerate, spaln or gth)
my $priority = 4;      # priority for hints
my $source   = "P";    # source for extrinsic file
my %genome;            # hash to store genome sequence
my $CDS; # output CDS instead of CDSpart hints
my $genome_file;
my $help;

if ( $#ARGV < 1 || $help ) {
    print $usage;
    exit;
}

GetOptions(
    'in=s'             => \$alignfile,
    'out=s'            => \$hintsfilename,
    'CDSpart_cutoff:i' => \$CDSpart_cutoff,
    'dir=s'            => \$dir,
    'maxintronlen:i'   => \$maxintronlen,
    'minintronlen:i'   => \$minintronlen,
    'priority:i'       => \$priority,
    'prg=s'            => \$prgsrc,
    'CDS!'             => \$CDS,
    'help!'            => \$help,
    'source:s'         => \$source,
    'genome_file=s'    => \$genome_file,
    'version!'         => \$printVersion
);

if ($printVersion) {
    print "align2hints.pl version $version\n";
    exit(0);
}

# mainly for usage within GALBA
if ( !defined($dir) ) {
    $dir = cwd();
}
my $last_char = substr( $dir, -1 );
if ( $last_char eq "\/" ) {
    chop($dir);
}

if ( defined($alignfile) ) {

    # check whether alignment file exists
    if ( !-e $alignfile ) {
        print STDERR
            "ERROR: Alignment file $alignfile does not exist. Please check.\n";
        exit(1);
    }
    else {
        $alignfile = rel2abs($alignfile);
    }
}

if ( !defined($prgsrc) ) {
    print STDERR
        "ERROR: Please assign the source programme with --prg. Possible "
        . "Options are 'gth' or 'miniprot'.\n";
    exit(1);
}

# check program source option
if (   $prgsrc ne "miniprot"
    && $prgsrc ne "gth")
{
    print STDERR
        "ERROR: Invalid value '$prgsrc' for option --prg. Possible Options "
        . "are 'miniprot' and 'gth'.\n";
    exit(1);
}

if ($prgsrc eq "miniprot") {
    $prgsrc = "miniprot2h";
}

if ( $prgsrc eq "gth" ) {
    $prgsrc           = "gth2h";
    $intron_threshold = $intron_threshold_gth;
}

if ( not( $prgsrc eq "miniprot2h" ) && defined($genome_file) )
{
    print STDERR
        "ERROR: program name is $prgsrc and a genome file was specified. "
        . "Will ignore genome file.\n";
} elsif ( $prgsrc eq "miniprot2h" && defined($genome_file) ) {
    open( GENOME, "<", $genome_file )
        or die("Could not open genome fasta file $genome_file!\n");
    my $header;
    while (<GENOME>) {
        chomp;
        if (m/^>(.*)/) {
            $genome{$1} = "";
            $header = $1;
        }
        else {
            $genome{$header} .= $_;
        }
    }
    close(GENOME) or die("Could not close genome fasta file $genome_file!\n");
}


if ($CDS) {
    $CDSpartid = "CDS";
}

open( ALN,   "<$alignfile" )     or die("Cannot open file: $alignfile\n");
open( HINTS, ">$hintsfilename" ) or die("Cannot open file: $hintsfilename");

while (<ALN>) {
    # skip if line contain star symbol
    if( m/\*/ ) {
        next;
    }
    chomp;
    my @f = split( /\t/, $_ );
    next unless ( scalar(@f) >= 8 );
    my $seqname = $f[0];
    my $type    = $f[2];
    my $start   = $f[3];
    my $end     = $f[4];
    my $score   = $f[5];
    my $strand  = $f[6];
    my $frame   = $f[7];
    if ( $end < $start ) {
        my $tmp = $start;
        $start = $end;
        $end   = $tmp;
    }

    # get target protein for gth
    if ( $type eq "mRNA" && $prgsrc eq "gth2h" ) {
        my @info   = split( /\=/, $f[8] );
        my @rnaid  = split( /;/,  $info[1] );    # $rnaid[0]
        my @geneid = split( /;/,  $info[2] );    # $geneid[0]
        @info = split( /\s/, $info[-1] );
        $parent
            = $info[0] . "_" . $seqname . "_" . $rnaid[0] . "_" . $geneid[0];
    }

    # define parent for miniprot
    if ($type eq "gene" && $prgsrc eq "miniprot2h"){
        $f[8] =~ m/gene_id \"([^"]+)\"/;
        $parent = $1;
    }

    # create start and stop hints from miniprot, stop codon excluded from CDS
    if ( ( $type eq "gene" && $prgsrc eq "miniprot2h" && defined($genome_file) ) )
    {
        my $pot_start;
        if ( $strand eq "+" ) {
            $pot_start = substr( $genome{$seqname}, $start - 1, 3 );
        }
        elsif ( $strand eq "-" ) {
            $pot_start = substr( $genome{$seqname}, $end - 3, 3 );
            $pot_start =~ tr/acgtACGT/tgcaTGCA/;
            $pot_start = reverse($pot_start);
        }
        if ( defined($pot_start) ) {
            if ( $pot_start =~ m/(ATG)|(TTG)|(GTG)|(CTG)/i ) {
                print_start( $seqname, $strand, $start, $end , $parent);
            }
        }
        my $pot_stop;
        if ( $strand eq "+" ) {
            $pot_stop = substr( $genome{$seqname}, $end, 3 );
        }
        else {
            $pot_stop = substr( $genome{$seqname}, $start - 4, 3 );
            $pot_stop =~ tr/acgtACGT/tgcaTGCA/;
            $pot_stop = reverse($pot_stop);
        }
        if ( $pot_stop =~ m/(TAA)|(TGA)|(TAG)/i ) {
            print_stop( $seqname, $strand, $start, $end );
        }
    }

    if ( $type eq "CDS" || $type eq "cds" || $type eq "protein_match" ) {
        if ( !$CDS ) {

            # CDSpart hint
            $start += $CDSpart_cutoff;
            $end -= $CDSpart_cutoff;
        }
        if ( $start > $end ) {
            $start = $end = int( ( $start + $end ) / 2 );
        }
        print HINTS
            "$seqname\t$prgsrc\t$CDSpartid\t$start\t$end\t$score\t$strand\t"
            . "$frame\tsrc=$source;grp=$parent;pri=$priority\n";
        if ( $prgsrc eq "miniprot2h" ) {
            get_intron( \@f );
        }
    }

    if ( $type eq "exon" && $prgsrc eq "gth2h" ) {
        get_intron( \@f );
    }
}
close(ALN)   or die("Could not close file $alignfile!\n");
close(HINTS) or die("Could not close file $hintsfilename!\n");

# intron hints for miniprot and gth
sub get_intron {
    my $line = shift;
    if($prgsrc eq "gth2h"){
        $intron_score = $prevScore + @{$line}[5] / 2;
    }else{
        $intron_score = 0;
    }
    if ( $prevParent ne $parent ) {
        if ( @{$line}[6] eq "-"
            && ( $prgsrc eq "miniprot2h" ) )
        {
            $intron_end = @{$line}[3] - 1
                ;   # these spliced aligners output in reverse order in genome
        }
        else {
            $intron_start = @{$line}[4] + 1;
        }
    }
    else {
        if ( @{$line}[6] eq "-"
            && ( $prgsrc eq "miniprot2h" ) )
        {
            $intron_start = @{$line}[4] + 1;
        }
        else {
            $intron_end = @{$line}[3] - 1;
        }
        if ( $intron_end < $intron_start ) {
            my $tmp = $intron_start;
            $intron_start = $intron_end;
            $intron_end   = $tmp;
        }

        if( $prgsrc eq "miniprot2h"){
            print HINTS
                "@{$line}[0]\t$prgsrc\tintron\t$intron_start\t$intron_end"
                . "\t$intron_score\t@{$line}[6]\t.\tsrc=$source;"
                . "grp=$parent;pri=$priority\n";
        }else{
            # for GTH check conditions: length of intron is at least $minintronlen and maximal
            # $maxintronlen and its score is greater than $intron_threshold
            if (   $intron_end - $intron_start + 1 >= $minintronlen
                && $intron_end - $intron_start + 1 <= $maxintronlen
                && ( !defined($intron_threshold)
                    || $intron_score > $intron_threshold )
                )
            {
                if ( ( defined($prevQend) && $prevQend + 1 == $qstart ) )
                {
                    print HINTS
                        "@{$line}[0]\t$prgsrc\tintron\t$intron_start\t$intron_end"
                        . "\t$intron_score\t@{$line}[6]\t.\tsrc=$source;"
                        . "grp=$parent;pri=$priority\n";
                }
            }
        }
        if ( @{$line}[6] eq "-"
            && ( $prgsrc eq "miniprot2h" ) )
        {
            $intron_end = @{$line}[3] - 1;
        }
        else {
            $intron_start = @{$line}[4] + 1;
        }
    }
    $prevScore  = @{$line}[5] / 2;
    $prevParent = $parent;
    $prevQend   = $qend if ( defined($qend) );
}

sub print_start {
    my $seqname = shift;
    my $strand  = shift;
    my $start   = shift;
    my $end     = shift;
    my $parent  = shift;
    print HINTS "$seqname\t$prgsrc\tstart\t";
    if ( $strand eq "+" ) {
        print HINTS "$start\t" . ( $start + 2 );
    }
    else {
        print HINTS ( $end - 2 ) . "\t$end";
    }
    print HINTS "\t.\t$strand\t0\tsrc=$source;grp=$parent;pri=$priority\n";
}

sub print_stop {
    my $seqname = shift;
    my $strand  = shift;
    my $start   = shift;
    my $end     = shift;
    print HINTS"$seqname\t$prgsrc\tstop\t";
    if ( $strand eq "+" ) {
        print HINTS ( $end - 2 ) . "\t$end";
    }
    else {
        print HINTS "$start\t" . ( $start + 2 );
    }
    print HINTS "\t.\t$strand\t0\tsrc=$source;grp=$parent;pri=$priority\n";
}
