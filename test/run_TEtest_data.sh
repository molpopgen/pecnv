#!/bin/bash

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

##DOWNLOAD THE DATA
if [ $WGET -eq 0 ] && [ $CURL -eq 0 ]
then
    echo "Error: neither wget nor curl were found in your \$PATH"
    exit 10;
fi

if [ ! -e TEtestData.tar ]
then
    if [ $WGET -eq 1 ]
    then
	wget http://devlaeminck.bio.uci.edu/Data/TEtestData.tar
    else
	if [ $CURL -eq 1 ]
	then
	    curl http://devlaeminck.bio.uci.edu/Data/TEtestData.tar -o TEtestData.tar
	fi
    fi
fi

##UNPACK THE ARCHIVE
if [ ! -d TEtestData ]
then
    tar xf TEtestData.tar
fi

##MAKE THE INPUT FILE
paste <(ls *_1.fastq.gz) <(ls *_2.fastq.gz) > INFILE

##Run the pecnv pipeline
pecnv.sh -i INFILE -r dmel-all-chromosome-r5.1_simplenames.fasta -c $CPU

