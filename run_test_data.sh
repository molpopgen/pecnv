#!sh

wget http://hpc.oit.uci.edu/~krthornt/references/dyak/dyak-all-chromosome-r1.3-newnames.fasta

bwa index dyak-all-chromosome-r1.3-newnames.fasta

wget http://hpc.oit.uci.edu:~/krthornt/dyak_genomic_data/NY4_*

ls -1 *.fastq.gz | sort > infile

pecnv.pl -infile infile -ref dyak-all-chromosome-r1.3-newnames.fasta -sample 0