#!/usr/bin/perl

use strict;
use Getopt::Long;

##Establish options

my $OUTDIR = "pecnv_output";
my $SAMPLEID = 0;
my $MINQUAL = 30;
my $MISMATCHES = 3;
my $GAPS = 0;
my $SAMPLES;
my $REFERENCE;
my $CPU=32;
GetOptions("outdir=s" => \$OUTDIR,
	   "minqual=i" => \$MINQUAL,
	   "mismatches=i" => \$MISMATCHES,
	   "gaps=i" => \$GAPS,
	   "infile=s" => \$SAMPLES,
	   "ref=s" => \$REFERENCE,
	   "sample=i" => \$SAMPLEID,
	   "cpu=i" => \$CPU
	  );

##Check that options are ok
if( ! defined( $SAMPLES ) )
  {
    die "Error, please specify input file with -infile <filename> option\n";
  }

if (! defined( $REFERENCE ) )
  {
    die "Error, please specify reference fasta file with -ref <filename>\n";
  }

##Make output dir

if ( ! -d $OUTDIR )
  {
    mkdir $OUTDIR;
  }

#Get FASTQ file names
my @FASTQFILES = ();
my @FASTQIDS=();

open I, "< $SAMPLES" or die "Error, could not open $SAMPLES for reading\n";

my $DUMMY = 0;
while(my $line = <I>)
  {
    chomp $line;
    push(@FASTQFILES,$line);
    push(@FASTQIDS,int($DUMMY/2));
    ++$DUMMY;
  }

##Rename the reads and align them
my $READ_DIR = 0;
for(my $i = 0 ; $i <= $#FASTQFILES ; ++$i )
  {
      system(qq{fastq_to_table $SAMPLEID $FASTQIDS[$i] $READ_DIR $FASTQFILES[$i] $OUTDIR/readfile.$i.fastq.gz});
      system(qq{bwa aln -t $CPU -l 13 -m 5000000 -I -R 5000 $REFERENCE $OUTDIR/readfile.$i.fastq.gz > $OUTDIR/readfile.$i.sai 2> $OUTDIR/alignment_stderr.$i});
    $READ_DIR = int(!$READ_DIR);
  }

##PE resolution step
for(my $i = 0 ; $i <= $#FASTQFILES ; $i += 2 )
  {
    my $sai1 = "$OUTDIR/readfile.$i.sai";
    my $j = $i+1;
    my $sai2 = "$OUTDIR/readfile.$j.sai";

    #Make position-sorted bamfile
    system(qq{bwa sampe -a 5000 -N 5000 -n 500 $REFERENCE $sai1 $sai2 $OUTDIR/readfile.$i.fastq.gz $OUTDIR/readfile.$j.fastq.gz 2> $OUTDIR/sampe_stderr.$i | samtools view -bS - 2> /dev/null | samtools sort -m 10000000 - $OUTDIR/bamfile.$FASTQIDS[$i]});
    unlink(qq{$sai1 $sai2});
  }

##Merge bam files
opendir( my $DH, $OUTDIR );
my @POSBAMS = grep { /\.bam$/ } readdir $DH;
closedir ($DH);

system(join(" ",qq{samtools merge $OUTDIR/merged_pos_sorted.bam },join(" ",@POSBAMS)));


#Make read name sorted bamfile.  We run this 2x b/c sometimes this step doesn't work
system(qq{samtools sort -n -m 10000000 $OUTDIR/merged_pos_sorted.bam $OUTDIR/temp});
system(qq{samtools sort -n -m 10000000 $OUTDIR/temp.bam $OUTDIR/merged_readsorted});
unlink(qq{$OUTDIR/temp.bam});

#Get distr. of mapping distances
system(qq{samtools view -f 2 $OUTDIR/merged_readsorted.bam | bwa_mapdistance $OUTDIR/mdist.gz});

#Get 99.9th quantileof mapping distances
system(qq{R --no-save --slave --args $OUTDIR/mdist.gz $OUTDIR/mquant.txt < mquant.R});

#Identify unusual read pairings
system(qq{samtools view -f 1 $OUTDIR/merged_readsorted.bam | bwa_bam_to_mapfiles $OUTDIR/cnv_mappings $OUTDIR/um});

#Finally, cluster the CNVs
system(qq{cluster_cnv $MINQUAL $MISMATCHES $GAPS $OUTDIR/div.gz  $OUTDIR/par.gz  $OUTDIR/ul.gz $OUTDIR/cnv_mappings_left.gz $OUTDIR/cnv_mappings_right.gz});
