#!/usr/bin/perl

use strict;
use Getopt::Long;
use Cwd;

my $WDIR = getcwd;

##Establish options and their default values
my $OUTDIR = "pecnv_output";
my $SAMPLEID = 0;
my $MINQUAL = 30;
my $MISMATCHES = 3;
my $GAPS = 0;
my $SAMPLES;
my $REFERENCE;
my $CPUMIN=8;
my $CPUMAX=64;
my $QUEUE;
my $JNAME = "pecnv";
GetOptions("outdir=s" => \$OUTDIR,
	   "minqual=i" => \$MINQUAL,
	   "mismatches=i" => \$MISMATCHES,
	   "gaps=i" => \$GAPS,
	   "infile=s" => \$SAMPLES,
	   "ref=s" => \$REFERENCE,
	   "sample=i" => \$SAMPLEID,
	   "cpumin=i" => \$CPUMIN,
	   "cpumin=i" => \$CPUMAX,
	   "q=s" => \$QUEUE,
	   "N=s" => \$JNAME
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

if (! defined( $QUEUE ) )
{
    die "Error, please specify queue with -q <queuename>\n";
}

##Make output dir
if ( ! -d $OUTDIR )
  {
    mkdir $OUTDIR;
  }

sub make_prefix
{
    if ( $_[0] < 10 )
    {
	return join("","0","$_[0]");
    }
    return $_[0];
}

##Count lines in infile
open I, "< $SAMPLES" or die "Error: could not open $SAMPLES for reading\n";

my $NFQ=0;

my @FASTQFILES=();
while(my $line =<I>)
{
    chomp $line;
    push(@FASTQFILES,$line);
    ++$NFQ;
}

close I;

my $SCRIPTNO = 0;

my $PREFIX = make_prefix($SCRIPTNO);

my $SNAME = join("",$PREFIX,"_rename_reads.sh");
my $JNAME_T = join("",$JNAME,"_rename");

##Job script to run fastq_to_table
open O, "> $SNAME";

my $HERE = <<HERE;
#!sh

#\$ -q $QUEUE
#\$ -t 1-$NFQ
#\$ -tc $CPUMAX
#\$ -N $JNAME_T
cd $WDIR

FILE=`head -n \$SGE_TASK_ID $SAMPLES | tail -n 1`

INDEX=\$((\$SGE_TASK_ID - 1))

READDIR=`expr \$INDEX % 2`

fastq_to_table $SAMPLEID \$INDEX \$READDIR \$FILE $OUTDIR/readfile.\$SGE_TASK_ID.fastq.gz

HERE

print O $HERE;

close O;

##Jobs to align
my $JNAME_O = $JNAME_T;

++$SCRIPTNO;
$PREFIX = make_prefix($SCRIPTNO);
$SNAME = join("",$PREFIX,"_align.sh");
$JNAME_T = join("",$JNAME,"_align");
open O, "> $SNAME";

$HERE=<<HERE;
#!sh

#\$ -q $QUEUE
#\$ -t 1-$NFQ
#\$ -pe openmp $CPUMIN-$CPUMAX
#\$ -N $JNAME_T
#\$ -hold_jid $JNAME_O;
cd $WDIR
bwa aln -t \$CORES -l 13 -m 5000000 -I -R 5000 $REFERENCE $OUTDIR/readfile.\$SGE_TASK_ID.fastq.gz > $OUTDIR/readfile.\$SGE_TASK_ID.sai 2> $OUTDIR/alignment_stderr.\$SGE_TASK_ID});
HERE
    
print O $HERE;

close O;

##PE resolution
++$SCRIPTNO;
$PREFIX = make_prefix($SCRIPTNO);
$SNAME = join("",$PREFIX,"_resolve.sh");
$JNAME_O=$JNAME_T;
$JNAME_T = join("",$JNAME,"_resolve");
open O, "> $SNAME";

open O, "> $SNAME";

my $NBAM = $NFQ/2;
$HERE=<<HERE;
#!sh

#\$ -q $QUEUE
#\$ -t 1-$NBAM
#\$ -pe openmp $CPUMIN-$CPUMAX
#\$ -N $JNAME_T
#\$ -hold_jid $JNAME_O;
cd $WDIR

INDEX=\$((\$SGE_TASK_ID-1)
FIRST=\$((2*(\$INDEX - 1)))
SECOND=\$((\$FIRST + 1))
bwa sampe -a 5000 -N 5000 -n 500 $REFERENCE $OUTDIR/readfile.\$FIRST.sai $OUTDIR/readfile.\$SECOND.sai $OUTDIR/readfile.\$FIRST.fastq.gz $OUTDIR/readfile.\$SECOND.fastq.gz 2> $OUTDIR/sampe_stderr.\$SGE_TASK_ID | samtools view -bS - 2> /dev/null | samtools sort -m 500000000 - $OUTDIR/bamfile.\$SGE_TASK_ID
rm -f $OUTDIR/readfile.\$FIRST.sai $OUTDIR/readfile.\$SECOND.sai
HERE
    
print O $HERE;

close O;




