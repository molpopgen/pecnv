#!/bin/sh

#$ -q krt,bio
#$ -pe openmp 64

cd /bio/krthornt/test_pecnv/pecnv

##Assumes bwa 0.5.9 and samtools are in your $PATH.
##You may want to "module load" if they are not. KRT keeps the versions used in Rogers et al. in his ~/bin
module load perl
module load boost/1.53.0
module load R

##Copy data from public-www as symlinks
#get unrenamed version of reference
cp -s /bio/public-www/krthornt/references/dyak/dyak-all-chromosome-r1.3-newnames.fasta .
#wget ftp://ftp.flybase.net/genomes/Drosophila_yakuba/dyak_r1.3_FB2011_08/fasta/dyak-all-chromosome-r1.3.fasta.gz
gunzip dyak-all-chromosome-r1.3.fasta.gz
bwa index dyak-all-chromosome-r1.3.fasta

cp -s /bio/public-www/krthornt/dyak_genomic_reads/NY42* .

ls -1 *.fastq.gz | sort > infile

#run pipeline
pecnv.pl -infile infile -ref dyak-all-chromosome-r1.3.fasta -sample 0

