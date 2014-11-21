#pecnv test scripts

This project contains two executable bash scripts.  One of these tests the CNV calling workflow using a subset of data from Rogers _et al._ (2014) PMID 24710518.  The other tests the transposable element (TE) calling pipeline from Cridlan et al. (2013) PMID 23883524.  The test of the TE workflow uses simulated data.

Please make sure that pecnv is installed on your system before running these tests.

The output of each script is several gigabytes of data!

##run_test_data.sh

This script performs CNV calling on two lanes of paired-end Illumina data from an inbred isofemale strain of _Drosophila yakuba_.  The script is very simple, and does the following:

1. Checks that you have wget and/or curl on your system.
2. Downloads a reference genome and 4 FASTQ files from Kevin Thornton's lab server.
3. Runs pecnv.sh on the data.

To run the test:

```
./run_test_data.sh
```

This script requires the bash shell, so please make sure that it is available.

##run_TE_test_data.sh

This script downloads a .tar archive from the Thornton lab server.  This archive contains a directory called TEtestData containing the following files:

```
TEtestData/dmel-all-chromosome-r5.1_simplenames.fasta
TEtestData/TE_position_r5.1
TEtestData/TE_frequencies_r5.1_chromnames
TEtestData/TE_frequencies_r5.1
TEtestData/line99_truth
TEtestData/line99_1.fastq.gz
TEtestData/line99_2.fastq.gz
TEtestData/transposon_all.fasta
TEtestData/TE_frequencies_r5.1_KRT
```

In order, these files are:

1. The _D. melanogaster_ r5.1 reference used in Cridland _et. al_ 2013.
2. The positions of annotated TEs in that reference genome.  The format is 3 columns: chromosomt, start, stop.  Start and stop count from position 1.
3. The contents of TE_position_r5.1 plus the TE label (Doc, roo, etc.) and the proportion of annotated TEs with that label.
4. An unsorted version of TE_position_r5.1.  You can ignore this file.
5. The insert site positions of 100 TEs that Julie Cridland randomly inserted into the reference genome
6. Lane 1 of data simulated from the reference genome
7. Lane 2 of data simulated from the reference genome
8. Sequences of Drosophila TEs in fasta format.  This is the FlyBase file + sequences for P and FP mined from Genbank. (The Dmel reference is P-free, so these sequences are included to aid annotation.)
9. The output of the one-liner shown below.  Same format as TE_position_r5.1

The TE_position_r5.1 file is a list of regions in the reference genome identified via blast searches to be similar to >= 75% of a sequence in transposon_all.fasta.  This is the file used in the Cridland _et al.__ work.

There are, however, lots of ways to make this file, such as RepeatMasker, etc.  One could also imagine other blastn-based criteria to take advantage of a known set of TE sequences.  For example, the following command will identify all regions in the reference genome hitting a sequence in transposon_all.fasta with an e-value of <= 1e-100 and then use BioConductor's [GenomicRanges](http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html) and [data.table](http://cran.r-project.org/web/packages/data.table/index.html)'s "fread" to greedily reduce the overlaps in the blast output and create a new file:

```
blastn -db dmel-all-chromosome-r5.1_simplenames.fasta -query transposon_all.fasta -outfmt 6 -num_threads 8 -evalue 1e-100  | cut -f 1,2,9,10 | awk '{if ($3 > $4) printf ("%s\t%s\t%s\t%s\n",$1,$2,$4,$3); else print $0}' | sort -k2,2 -k3,3n | Rscript -e 'library(GenomicRanges);library(data.table);x=fread("cat /dev/stdin");gr=GRanges(seqnames=x$V2,ranges=IRanges(start=x$V3,end=x$V4));grr=as.data.frame(reduce(gr));write.table(cbind(as.character(grr$seqnames),grr$start,grr$end),file="TE_position_r5.1_KRT",row.names=F,col.names=F,quote=F)'
```

Let's decode that one-liner, as it is tough to read.  The steps are:

1. Run blastn of the TE fasta file against the reference genome, printing tabular output to STDOUT
2. Extract query name, hit name, and start/stop positions along hit
3. For alignments on the minus strand, swap start stop if needed
4. Sort by chromosome name and position along chromosome
5. Load the data into R, create GenomicRanges, and then "reduce", which greedily collapses the overlaps.
6. Write the reduced data.table to a file called TE_position_r5.1_KRT

To run the script:

```
./run_TEtest_data.sh
```

What the script does:

1. Download and untar the archive
2. Create the infiles
3. Run pecnv.sh
4. Run teclust on the output of the pecnv.sh pipeline using the Cridland et al. method
5. If phrap and GNU parallel are available on your system, it'll perform _de novo_ assembly of lots of events

Once you have run this script, you may be interested in running TEvariationsOnTheme.sh, which illustrates different ways to run teclust.