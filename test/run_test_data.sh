#!/usr/bin/env bash

CPU=4

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

if [ $WGET -eq 0 ] && [ $CURL -eq 0 ]
then
    echo "Error: neither wget nor curl were found in your \$PATH"
    exit 10;
fi

if [ $WGET -eq 1 ]
then
    if [ ! -e dyak-all-chromosome-r1.3-newnames.fasta ]
	then
	wget http://devlaeminck.bio.uci.edu/Data/references/dyak-all-chromosome-r1.3-newnames.fasta
    fi
    for READFILE in NY42_06_21_2010_54_1.fastq.gz NY42_06_21_2010_54_2.fastq.gz NY42_09_07_2010_75_1.fastq.gz NY42_09_07_2010_75_2.fastq.gz
    do
	if ! [ -e $READFILE ]
	then
	    wget http://devlaeminck.bio.uci.edu/Data/dyak_genomic_reads/$READFILE
	fi
    done
else
    if [ ! -e dyak-all-chromosome-r1.3-newnames.fasta ]
    then
	curl http://devlaeminck.bio.uci.edu/references/dyak-all-chromosome-r1.3-newnames.fasta -o dyak-all-chromosome-r1.3-newnames.fasta
    fi
#Download two lanes of data for 1 sample--that's enough to see how things progress
    for READFILE in NY42_06_21_2010_54_1.fastq.gz NY42_06_21_2010_54_2.fastq.gz NY42_09_07_2010_75_1.fastq.gz NY42_09_07_2010_75_2.fastq.gz
    do
	if ! [ -e $READFILE ]
	then
	    curl http://devlaeminck.bio.uci.edu/Data/dyak_genomic_reads/$READFILE -o $READFILE
	fi
    done
fi

#bwa index dyak-all-chromosome-r1.3-newnames.fasta

#Make the input file. 
#This command is specific to these data
#The format is two columns, separated by a SINGLE TAB
#LEFT READ FILE (tab) RIGHT READ FILE
#etc. for all lanes
paste <(ls *_1.fastq.gz) <(ls *_2.fastq.gz) > INFILE

pecnv.sh -i INFILE -r dyak-all-chromosome-r1.3-newnames.fasta -c $CPU
