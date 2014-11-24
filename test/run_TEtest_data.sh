#!/usr/bin/env bash

CPU=4
PHRAPCPU=8

command_exists () {
    type "$1" &> /dev/null ;
}

WGET=0
CURL=0

if command_exists wget
then
    WGET=1
fi

if command_exists curl
then 
    CURL=1
fi

##DOWNLOAD THE DATA
if [ $WGET -eq 0 ] && [ $CURL -eq 0 ]
then
    echo "Error: neither wget nor curl were found in your \$PATH"
    exit 10;
fi

TARFILE=TEtestData.0.1.5.tar

if [ ! -e $TARFILE ]
then
    if [ $WGET -eq 1 ]
    then
	wget http://devlaeminck.bio.uci.edu/Data/$TARFILE
    else
	if [ $CURL -eq 1 ]
	then
	    curl http://devlaeminck.bio.uci.edu/Data/$TARFILE -o $TARFILE
	fi
    fi
fi

##UNPACK THE ARCHIVE
if [ ! -d TEtestData ]
then
    tar xf $TARFILE
fi

cd TEtestData

##MAKE THE INPUT FILE
paste <(ls *_1.fastq.gz) <(ls *_2.fastq.gz) > INFILE

##Run the pecnv pipeline
pecnv.sh -S line99 -i INFILE -r dmel-all-chromosome-r5.1_simplenames.fasta -c $CPU

##Get the 99th quantile of the insert size distribution
u99=`pecnv_insert_qtile pecnv_output/pecnv_bamfile.mdist.gz 0.99`

##Run teclust
teclust -s line99 -b pecnv_output/pecnv_bamfile_sorted.bam -t TE_position_r5.1.bed -o teclust_output.gz -u pecnv_output/pecnv_bamfile.um_u.csv.gz -m pecnv_output/pecnv_bamfile.um_m.csv.gz -i $u99 -p phrapdir

##If we have the right stuff on the system, assemble in parallel using phrap
CANASSEMBLE=1
for needed in parallel phrap
do
    PM=`which $needed`
    if [ -z ${PM} ]
    then
	>&2 echo "Cannot find $needed in your user's PATH"
	CANASSEMBLE=0
    else
	echo "$needed is $PM"
    fi
done

if [ $CANASSEMBLE -eq 1 ]
then
    find phrapdir -name "*.fasta" | parallel --jobs $PHRAPCPU --timeout 300 "phrap {} -vector_bound 0 -forcelevel 10 -minscore 10 -minmatch 10 -new_ace 2> {}.stderr > {}.stdout"
fi
