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
    >&2 echo "Example:"
    >&2 echo "$0 -i readfile.txt -r reference.fa"
    exit 1
}

OUTDIR="pecnv_output";
SAMPLEID=0;
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
	-- ) shift; break ;;
    * ) break ;;
  esac
done

##These params are important for TE detection (Cridland et al.)
BWAEXTRAPARMS="-l 13 -m $ALNMEM -I -R 5000"

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

##VALIDATE THE INPUT PARAMS
if [ -z ${SAMPLES+x} ]; then >&2 echo "Error: no input file specified"; usage; else echo "Input file name is set to '$SAMPLES'"; fi
if [ -z ${REFERENCE+x} ]; then >&2 echo "Error: no reference file specified"; usage; else echo "Reference file name is set to '$REFERENCE'"; fi

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

NPAIRS=${#LEFTS[@]}
PAIR=0
while [ $PAIR -lt $NPAIRS ]
do
    LEFTREADS=${LEFTS[$PAIR]}
    RIGHTREADS=${RIGHTS[$PAIR]}
    BAMFILEBASE="$OUTDIR/intermediate_bamfile$PAIR"
    if [ $NPAIRS -eq 1 ]
    then
	BAMFILEBASE="$OUTDIR/$BAMFILESTUB"
    fi
    #ALIGN THIS PAIR 
    #We go straight to sorted BAM output via process substitution
    bwa sampe $REFERENCE <(bwa aln -t $CPU $BWAEXTRAPARMS $REFERENCE $LEFTREADS 2> $OUTDIR/$align_stderr_1.$PAIR) <(bwa aln -t $CPU $BWAEXTRAPARMS $REFERENCE $RIGHTREADS 2> $OUTDIR/$align_stderr_1.$PAIR) $LEFTREADS $RIGHTREADS 2> $OUTDIR/sampe.$PAIR.stderr | samtools view -bS - 2> $OUTDIR/view.$PAIR.stderr | samtools sort -m $SORTMEM - $BAMFILEBASE 2> $OUTDIR/sort.$PAIR.stderr
    PAIR=$(($PAIR+1))
done

###2. Merge bam files if necessary
if [ $NPAIRS -gt 1 ]
then
    samtools merge $OUTDIR/"$BAMFILESTUB"_sorted.bam $OUTDIR/intermediate_bamfile*.bam
fi

###3. Sort bam file by read name.  This is done 2x b/c samtools sorting on read name has been unreliable in the past
samtools sort -n -m $SORTMEM $OUTDIR/"$BAMFILEBASE"_sorted.bam $OUTDIR/temp
samtools sort -n -m $SORTMEM $OUTDIR/temp.bam $OUTDIR/"$BAMFILESTUB"_readsorted

if [ -s $OUTDIR/"$BAMFILESTUB"_readsorted ]
then
    rm -f $OUTDIR/temp*.bam
fi


###4. Collect unusual read pairings
samtools view -f 1 $OUTDIR/"$BAMFILESTUB"_readsorted.bam | process_readmappings $OUTDIR/cnv_mappings $OUTDIR/um $OUTDIR/mdist.gz

###5. Get quantile of mapping distance
Rscript -e "x=read.table(\"$OUTDIR/mdist.gz\",header=T);z=which(x\$cprob >= 0.999);y=x\$distance[z[1]];write(y,\"$OUTDIR/mquant.txt\")"

MD=`head -n 1 $OUTDIR/mquant.txt`

###6. Cluster
cluster_cnv $MINQUAL $MISMATCHES $GAPS $MD $OUTDIR/div_clusters.gz  $OUTDIR/par_clusters.gz  $OUTDIR/ul_clusters.gz $OUTDIR/cnv_mappings.csv.gz
