#!sh

#Change this if your system has < 32 CPU
CPU=32

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
    wget http://hpc.oit.uci.edu/~krthornt/references/dyak/dyak-all-chromosome-r1.3-newnames.fasta
else
    curl http://hpc.oit.uci.edu/~krthornt/references/dyak/dyak-all-chromosome-r1.3-newnames.fasta -o dyak-all-chromosome-r1.3-newnames.fasta
fi

bwa index dyak-all-chromosome-r1.3-newnames.fasta

#Download two lanes of data for 1 sample--that's enough to see how things progress
wget http://hpc.oit.uci.edu/~krthornt/dyak_genomic_reads/NY42_06_21_2010_54_1.fastq.gz
wget http://hpc.oit.uci.edu/~krthornt/dyak_genomic_reads/NY42_06_21_2010_54_2.fastq.gz
wget http://hpc.oit.uci.edu/~krthornt/dyak_genomic_reads/NY42_09_07_2010_75_1.fastq.gz
wget http://hpc.oit.uci.edu/~krthornt/dyak_genomic_reads/NY42_09_07_2010_75_2.fastq.gz

ls -1 *.fastq.gz | sort > infile

pecnv.pl -infile infile -ref dyak-all-chromosome-r1.3-newnames.fasta -sample 0 -cpu $CPU