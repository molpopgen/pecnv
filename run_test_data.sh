#!sh

#Change this if your system has < 32 CPU
CPU=32

wget http://hpc.oit.uci.edu/~krthornt/references/dyak/dyak-all-chromosome-r1.3-newnames.fasta

bwa index dyak-all-chromosome-r1.3-newnames.fasta

#Download two lanes of data for 1 sample--that's enough to see how things progress
wget http://hpc.oit.uci.edu/~krthornt/dyak_genomic_reads/NY42_06_21_2010_54_1.fastq.gz
wget http://hpc.oit.uci.edu/~krthornt/dyak_genomic_reads/NY42_06_21_2010_54_2.fastq.gz
wget http://hpc.oit.uci.edu/~krthornt/dyak_genomic_reads/NY42_09_07_2010_75_1.fastq.gz
wget http://hpc.oit.uci.edu/~krthornt/dyak_genomic_reads/NY42_09_07_2010_75_2.fastq.gz

ls -1 *.fastq.gz | sort > infile

pecnv.pl -infile infile -ref dyak-all-chromosome-r1.3-newnames.fasta -sample 0 -cpu $CPU