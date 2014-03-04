#!/usr/bin/perl

use strict;
use Getopt::Long;
use Parallel::ForkManager;

##Establish options and their default values
my $OUTDIR = "pecnv_output";
my $SAMPLEID = 0;
my $MINQUAL = 30;
my $MISMATCHES = 3;
my $GAPS = 0;
my $SAMPLES;
my $REFERENCE;
my $CPU=32;
my $HAVEROCKSORT;
my $ROCKMEM="32G";
GetOptions("outdir=s" => \$OUTDIR,
	   "minqual=i" => \$MINQUAL,
	   "mismatches=i" => \$MISMATCHES,
	   "gaps=i" => \$GAPS,
	   "infile=s" => \$SAMPLES,
	   "ref=s" => \$REFERENCE,
	   "sample=i" => \$SAMPLEID,
	   "cpu=i" => \$CPU,
	   "rocksort" => \$HAVEROCKSORT,
	   "m=s" => \$ROCKMEM
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

my $ROCKMEMVAL=$ROCKMEM;
my $ROCKMEMUNIT=$ROCKMEM;
$ROCKMEMVAL =~ s/[A-Za-z]+//go;
$ROCKMEMUNIT =~ =s/[0-9]+//go;

my $ROCKMEMCPU = int($ROCKMEMVAL/$CPU);


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

##Rename the reads
#my $READ_DIR = 0;
#my @CONVERTCLI=();
#for(my $i = 0 ; $i <= $#FASTQFILES ; ++$i )
#  {
#      push(@CONVERTCLI,qq{fastq_to_table $SAMPLEID $FASTQIDS[$i] $READ_DIR $FASTQFILES[$i] $OUTDIR/readfile.$i.csv.gz $OUTDIR/readfile.$i.fastq.gz});
#      $READ_DIR = int(!$READ_DIR);
#  }

my $pm = new Parallel::ForkManager($CPU);

#foreach my $C (@CONVERTCLI)
#{
#    $pm->start and next;
#    system(qq{$C});
#    $pm->finish;
#}

#$pm->wait_all_children;

##Align the reads
for(my $i = 0 ; $i <= $#FASTQFILES ; ++$i )
  {
      system(qq{bwa aln -t $CPU -l 13 -m 5000000 -I -R 5000 $REFERENCE $FASTQFILES[$i] > $OUTDIR/readfile.$i.sai 2> $OUTDIR/alignment_stderr.$i});
#      $READ_DIR = int(!$READ_DIR);
  }

##PE resolution step
my @RESOLVE = ();
my @SAI=();
for(my $i = 0 ; $i <= $#FASTQFILES ; $i += 2 )
  {
    my $sai1 = "$OUTDIR/readfile.$i.sai";
    my $j = $i+1;
    my $sai2 = "$OUTDIR/readfile.$j.sai";
    push(@SAI,$sai1);
    push(@SAI,$sai2);
    #Make position-sorted bamfile
    if ( defined ($HAVEROCKSORT) )
    {
	push(@RESOLVE,qq{bwa sampe -a 5000 -N 5000 -n 500 $REFERENCE $sai1 $sai2 $FASTQFILES[$i] $FASTQFILES[$j] 2> $OUTDIR/sampe_stderr.$i | samtools view -bS - 2> /dev/null | samtools rocksort -m $ROCKMEM - $OUTDIR/bamfile.$FASTQIDS[$i]});
    }
    else
    {
	push(@RESOLVE,qq{bwa sampe -a 5000 -N 5000 -n 500 $REFERENCE $sai1 $sai2 $FASTQFILES[$i] $FASTQFILES[$j] 2> $OUTDIR/sampe_stderr.$i | samtools view -bS - 2> /dev/null | samtools sort -m 500000000 - $OUTDIR/bamfile.$FASTQIDS[$i]});
    }
  }

$pm->set_max_procs(int($CPU/2));

foreach my $R (@RESOLVE)
{
    $pm->start and next;
    system(qq{$R});
    $pm->finish;
}

$pm->wait_all_children;

##Cleanup the .sai files
foreach my $S (@SAI)
{
    unlink(qq{$S});
}

##Merge bam files
opendir( my $DH, $OUTDIR );
my @POSBAMS = grep { /\.bam$/ } readdir $DH;
closedir ($DH);

##Prepend output dir path to bamfile names
foreach my $P (@POSBAMS)
{
    $P = join("/",$OUTDIR,$P);
}

if ($#POSBAMS > 0)
{
    system(join(" ",qq{samtools merge $OUTDIR/merged_pos_sorted.bam },join(" ",@POSBAMS)));
}
else
{
    system(qq{mv $OUTDIR/bamfile.0.bam $OUTDIR/merged_pos_sorted.bam});
}

##Index.  This pipeline doesn't use the index, but it is good to have
system(qq{samtools index $OUTDIR/merged_pos_sorted.bam});

##Clean up
unlink(join(" ",@POSBAMS));

#Make read name sorted bamfile.  
if ( defined($HAVEROCKSORT) )
{
    system(qq{samtools rocksort -n -m $ROCKMEMCPU$ROCKMEMUNIT -@ $CPU $OUTDIR/merged_pos_sorted.bam $OUTDIR/merged_readsorted});
}
else
{
#We run this 2x b/c sometimes this step doesn't work and results in an icompletely-sorted bam file
    system(qq{samtools sort -n -m 500000000 $OUTDIR/merged_pos_sorted.bam $OUTDIR/temp});
    system(qq{samtools sort -n -m 500000000 $OUTDIR/temp.bam $OUTDIR/merged_readsorted});
    unlink(qq{$OUTDIR/temp.bam});
}

#Get distr. of mapping distances
system(qq{samtools view -f 2 $OUTDIR/merged_readsorted.bam | bwa_mapdistance $OUTDIR/mdist.gz});

#Get 99.9th quantile of mapping distances
#Executed via a HERE document so that we don't need path to an R script
system(qq{R --no-save --slave --args $OUTDIR/mdist.gz $OUTDIR/mquant.txt <<TEST;
n=commandArgs(trailing=TRUE)
x=read.table(n[1],header=T)
z=which(x\\\$cprob >= 0.999)
y=x\\\$distance[z[1]]
write(y,n[2])
TEST
});

#Identify unusual read pairings
system(qq{samtools view -f 1 $OUTDIR/merged_readsorted.bam | bwa_bam_to_mapfiles $OUTDIR/cnv_mappings $OUTDIR/um});

#Finally, cluster the CNVs
my @MDIST=`cat $OUTDIR/mquant.txt`;
my $MD = shift(@MDIST);
chomp $MD;

system(qq{cluster_cnv $MINQUAL $MISMATCHES $GAPS $MD $OUTDIR/div.gz  $OUTDIR/par.gz  $OUTDIR/ul.gz $OUTDIR/cnv_mappings_left.csv.gz});# $OUTDIR/cnv_mappings_right.csv.gz});
