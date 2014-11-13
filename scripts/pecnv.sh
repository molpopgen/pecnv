#!/bin/bash

usage(){
    >&2 echo "Usage: $0 options"
    >&2 echo "Mandatory options are:"
    >&2 echo " -i/--infile = A text file with the pairs of fastq file names listed.  One pair per line.  Must be TAB-separated"
    >&2 echo " -r/--reference = The reference file name to use for alignment."
    >&2 echo "Optional arguments are:"
    >&2 echo " -o/--outdir = output directory name for the output.  Will be created if it does not exist."
    >&2 echo " -q/--minqual = Min. mapping quality for clustering.  Default = 30"
    >&2 echo " -m/--mismatches = Max. number alignment mismatches to allow for clustering.  Default = 3"
    >&2 echo " -g/--gaps = Max. number of alignment gaps to allow for clustering.  Default = 0"
    >&2 echo " -c/--cpu = Number of CPU to use for BWA alignment.  This should be 1/2 the number of cores you have available,"
    >&2 echo "            because we align both reads from a lane at once, skipping intermediate files via process subsitution."
    >&2 echo "            Default = 32"
    >&2 echo " -s/--sortmem = Memory to be used by samtools sort.  Default = 50000000"
    >&2 echo " -a/--alnmem = Memory to be use by bwa aln step.  Default = 5000000"
    >&2 echo " -b/--bamfilebase = Prefix for bam file.  Default = pecnv_bamfile"
    >&2 echo " -u/--ulimit = max RAM usage for processing BAM file.  Unit is in gigabytes, e.g., 5 will be converted to 5*1024^2 bytes"
    >&2 echo "Example:"
    >&2 echo "$0 -i readfile.txt -r reference.fa"
    exit 1
}

OUTDIR="pecnv_output";
MINQUAL=30;
MISMATCHES=3;
GAPS=0;
CPU=32;
SORTMEM=50000000
ALNMEM=5000000
BAMFILESTUB="pecnv_bamfile"

while true; do
    case "$1" in
	-o | --outdir ) OUTDIR="$2"; shift 2 ;;
	-q | --minqual ) MINQUAL="$2"; shift 2;;
	-m | --mismatches ) MISMATCHES="$2"; shift 2;;
	-g | --gaps ) GAPS="$2"; shift 2;;
	-i | --infile ) SAMPLES="$2"; shift 2;;
	-c | --cpu ) CPU="$2"; shift 2;;
	-r | --reference ) REFERENCE="$2"; shift 2;;
	-s | --sortmem ) SORTMEM="$2"; shift 2;;
	-a | --alnmem ) ALNMEM="$2"; shift 2;;
	-b | --bamfilebase ) BAMFILESTUB="$2" ; shift 2;;
	-u | --ulimit ) MAXRAM=`echo "$2*1025^2"|bc -l` ; shift 2;;
	-- ) shift; break ;;
    * ) break ;;
  esac
done

##These params are important for TE detection (Cridland et al.)
BWAEXTRAPARMS="-l 13 -m $ALNMEM -I -R 5000"

##VALIDATE THE INPUT PARAMS
if [ -z ${SAMPLES+x} ]; then >&2 echo "Error: no input file specified"; usage; else echo "Input file name is set to '$SAMPLES'"; fi
if [ -z ${REFERENCE+x} ]; then >&2 echo "Error: no reference file specified"; usage; else echo "Reference file name is set to '$REFERENCE'"; fi

#Check for executable dependencies
for needed in bwa samtools process_readmappings Rscript cluster_cnv
do
    PM=`which $needed`
    if [ -z ${PM} ]
    then 
	>&2 echo "Error: $needed not found in your \$PATH"; 
	exit; 
    else 
	echo "$needed is $PM"; 
    fi;
done

if [ ! -e $SAMPLES ]
then
    >&2 echo "Error: file $SAMPLES does not exist as a regular file"
    exit
fi

if [ ! -e $REFERENCE ]
then
    >&2 echo "Error: file $REFERENCE does not exist as a regular file"
    exit
fi

if [ ! -e $REFERENCE.bwt ]
then
    >&2 echo "Reference does not seem to have been indexed.  Running bwa index..."
    bwa index $REFERENCE
    if [ ! -e $REFERENCE.bwt ]
    then
	>&2 echo "Error: I tried to index your reference, but still cannot find $REFERENCE.bwt.  Perhaps there is a permissions issue?"
	exit;
    fi
fi

if [ ! -d $OUTDIR ]
then
    echo "Making output directory $OUTDIR"
    mkdir $OUTDIR
fi

##THE WORK STARTS HERE

###1. Align read pairs and make bam files
LEFTS=(`cut -f 1 $SAMPLES`)
RIGHTS=(`cut -f 2 $SAMPLES`)

if [ ${#LEFTS[@]} -ne ${#RIGHTS[@]} ]
then
    >&2 echo "Error: unequal number of left and right read files in $SAMPLES"
    exit
fi

if [ ! -e $OUTDIR/"$BAMFILESTUB"_sorted.bam ]
then
    NPAIRS=${#LEFTS[@]}
    PAIR=0
    while [ $PAIR -lt $NPAIRS ]
    do
	LEFTREADS=${LEFTS[$PAIR]}
	RIGHTREADS=${RIGHTS[$PAIR]}
	BAMFILEBASE="$OUTDIR/intermediate_bamfile$PAIR"
	if [ $NPAIRS -eq 1 ]
	then
	    BAMFILEBASE="$OUTDIR/$BAMFILESTUB"_sorted
	fi
    #ALIGN THIS PAIR 
    #We go straight to sorted BAM output via process substitution
	bwa sampe -a 5000 -N 5000 -n 500 $REFERENCE <(bwa aln -t $CPU $BWAEXTRAPARMS $REFERENCE $LEFTREADS 2> $OUTDIR/align_stderr_1.$PAIR) <(bwa aln -t $CPU $BWAEXTRAPARMS $REFERENCE $RIGHTREADS 2> $OUTDIR/align_stderr_1.$PAIR) $LEFTREADS $RIGHTREADS 2> $OUTDIR/sampe.$PAIR.stderr | samtools view -bS - 2> $OUTDIR/view.$PAIR.stderr | samtools sort -m $SORTMEM - $BAMFILEBASE 2> $OUTDIR/sort.$PAIR.stderr
	PAIR=$(($PAIR+1))
    done

###2. Merge bam files if necessary
    if [ $NPAIRS -gt 1 ]
    then
	samtools merge $OUTDIR/"$BAMFILESTUB"_sorted.bam $OUTDIR/intermediate_bamfile*.bam
	rm -f $OUTDIR/intermediate_bamfile*.bam
    fi
else
    echo $OUTDIR/"$BAMFILESTUB"_sorted.bam exists, so skipping alignment step
fi

###3. Collect unusual read pairings and estimate insert size distributions
if [ -z ${MAXRAM+x} ]
then
    >&2 echo "Using default ulimit"
else
    ulimit -v $MM
fi

process_readmappings $OUTDIR/"$BAMFILESTUB"_sorted.bam $OUTDIR/$BAMFILESTUB.cnv_mappings $OUTDIR/$BAMFILESTUB.um 
bwa_mapdistance $OUTDIR/"$BAMFILESTUB"_sorted.bam $OUTDIR/$BAMFILESTUB.mdist.gz

###4. Cluster (uses Rscript to get the 99.9th quantile of insert size distribution)
cluster_cnv $MINQUAL $MISMATCHES $GAPS `pecnv_insert_qtile $OUTDIR/$BAMFILESTUB.mdist.gz 0.999` $OUTDIR/$BAMFILESTUB.div.gz  $OUTDIR/$BAMFILESTUB.par.gz  $OUTDIR/$BAMFILESTUB.ul.gz $OUTDIR/$BAMFILESTUB.cnv_mappings.csv.gz
