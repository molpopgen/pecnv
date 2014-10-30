#!/bin/bash

#Change this if your system has < 32 CPU
CPU=8

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
    wget http://devlaeminck.bio.uci.edu/Data/references/dyak-all-chromosome-r1.3-newnames.fasta
#Download two lanes of data for 1 sample--that's enough to see how things progress
    wget http://devlaeminck.bio.uci.edu/Data/dyak_genomic_reads/NY42_06_21_2010_54_1.fastq.gz
    wget http://devlaeminck.bio.uci.edu/Data/dyak_genomic_reads/NY42_06_21_2010_54_2.fastq.gz
    wget http://devlaeminck.bio.uci.edu/Data/dyak_genomic_reads/NY42_09_07_2010_75_1.fastq.gz
    wget http://devlaeminck.bio.uci.edu/Data/dyak_genomic_reads/NY42_09_07_2010_75_2.fastq.gz
else
    curl http://devlaeminck.bio.uci.edu/references/dyak-all-chromosome-r1.3-newnames.fasta -o dyak-all-chromosome-r1.3-newnames.fasta
#Download two lanes of data for 1 sample--that's enough to see how things progress
    curl http://devlaeminck.bio.uci.edu/Data/dyak_genomic_reads/NY42_06_21_2010_54_1.fastq.gz -o NY42_06_21_2010_54_1.fastq.gz
    curl http://devlaeminck.bio.uci.edu/Data/dyak_genomic_reads/NY42_06_21_2010_54_2.fastq.gz -o NY42_06_21_2010_54_2.fastq.gz
    curl http://devlaeminck.bio.uci.edu/Data/dyak_genomic_reads/NY42_09_07_2010_75_1.fastq.gz -o NY42_09_07_2010_75_1.fastq.gz
    curl http://devlaeminck.bio.uci.edu/Data/dyak_genomic_reads/NY42_09_07_2010_75_2.fastq.gz -o NY42_09_07_2010_75_2.fastq.gz
fi

bwa index dyak-all-chromosome-r1.3-newnames.fasta

#Make the input file. 
#This command is specific to these data
#The format is two columns, separated by a SINGLE TAB
#LEFT READ FILE (tab) RIGHT READ FILE
#etc. for all lanes
paste <(ls *_1.fastq.gz) <(ls *_2.fastq.gz) > INFILE

#ls -1 *.fastq.gz | sort > infile

#pecnv.pl -infile infile -ref dyak-all-chromosome-r1.3-newnames.fasta -sample 0 -cpu $CPU
pecnv.sh -i infile -r dyak-all-chromosome-r1.3-newnames.fasta -cpu $CPU