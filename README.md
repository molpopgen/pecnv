#Introduction
_pecnv_ is a collection of C++ programs used for the detection of copy-number variation (CNV) by means of paired-end mapping of short reads.  

This pipeline corresponds to the clustering algorithm described and used in

*  Rogers, R. L.,J. M. Cridland, L. Shao, T. T. Hu, P. Andolfatto and K. R. Thornton (2014) Landscape of standing variation for tandem duplications in _Drosophila yakuba_ and _Drosophila simulans_.  Molecular Biology and Evolution __31__: 1750-1766 PMID 24710518, [Manuscript](http://mbe.oxfordjournals.org/content/31/7/1750.abstract.html)

This software also a C++ implementation of the transposable element (TE) detection pipeline descrbed in:

* Cridland, J.M., S.J. MacDonald, A.D. Long, and K.R Thornton (2013) Abundance and Distribution of Transposable Elements in Two _Drosophila_ QTL Mapping Resources  Molecular Biology and Evolution __30__: 2311-2327. PMID 23883524 [Manuscript](http://mbe.oxfordjournals.org/content/30/10/2311.full)

The full version of the TE detection pipeline, along with test data, is available from the Thornton lab [website](http://www.molpopgen.org/tepipeline/line99_example.tar.gz).  Users interested primarily in TE detection are encouraged to download that archive and study how it works.

This code has been used by the lab for detecting tandem duplications from short-read data, and the accompanying perl/shell scripts are intended to walk a user through doing such an analysis.

#A note on version numbers

1.  Release 0.1.0 = the precise version used in Rogers et al. and Cridland et al.
2.  Release 0.1.1 = a bug fix resolving issue number 1 regarding inconsistent SAM flags within read pairs.  This bug does not affect the conclusions of the previous work because it affects very few reads that make it into the final clustering steps.
3.  Release 0.1.2 = No more need to rename references or read names.  This shaves a load of time off of the workflow.
4.  Release 0.1.3 = re-implementation of process_readmappings to read directly from BAM files.  No more need to make a BAM file sorted by read name.  This shaves hours off of the pipeline.

##Tentative plan for future releases

1.  0.1.4 will integrate the workflow for calling transposable-element insertions from PE mapping data using the Cridland _et al._ approach.  This version is likely to add a dependency on the [boost](http://www.boost.org) program_options library, in order to make simplify the interface to the TE calling programs.

#Installation

##Compile-time dependencies

This software requires the following libraries to compile:

1. [libsequence](http://www.github.com/molpopgen/libsequence) - version 1.8.3 or greater
2. [zlib](http://zlib.net) (This is a also a libsequence dependency)
3. [htslib](http://htslib.org)
4. [boost](http://www.boost.org) is used for the processing of command-line options for the TE detection program.  Specifically, you need the boost_program_options and boost_system libraries installed.

__Please make sure that your installation of libsequence makes use of htslib.  That requires that htslib be installed prior to installing libsequence, so that the latter will compile features depending on the former.__

##Run-time dependencies

In order to run the pecnv.sh script, you need the following tools on your system:

* The [bash](http://www.gnu.org/software/bash/) shell.  The pecnv.sh script is written in bash.
* The [Rscript](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html) interface to [R](http://www.r-project.org)

For the detection of putative transposable element insertions, the following programs are optional:

* [phrap](http://www.phrap.org), version >= 1.09 (it must be compatible with short-read data)
* [GNU parallel](http://www.gnu.org/software/parallel/)
* [NCBI blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi)

##Compilation and installation

For systems where all dependencies are "where they're supposed to be":

```
./configure
make
sudo make install
```

For systems where dependencies are in a custom location.  For this example, I assume that they are installed in your users' home folder:

```
./configure CXXFLAGS=-I$HOME/include LDFLAGS=-L$HOME/lib
make
sudo make install
```

To install into your user's home folder:

```
./configure --prefix=$HOME
make
make install
```

You may mix and match the above.

##What this package installs

This package installs several executables.  The two that users are likely to run directly are:

*pecnv.sh is an executable bash script that runs the CNV calling pipeline on paired-end short-read Illumina data (see next section)
*teclust is a C++ program implementing the TE detection method of Cridland _et al._ (2013) PMID 23883524.  This program is explained in more detail below.

The following executables are installed and called by pecnv.sh:
*process_readmappings is a C++ program that collects read pairs in unxepected orientations (divergent, parallel, unlinked, and unique/multi) from a BAM file
*bwa_mapdistance is a C++ program that summarizes the ECDF of the insert size distribution based on uniquely-mapped paired reads in a BAM file
*pecnv_insert_qtile is an executable Rscript that processes the output of bwa_mapdistance to estimate the upper quantile of the insert size distribtution
*cluster_cnv is a C++ program that clusters the output of process_readpmappings into putative CNV calls

"Power users" are encouraged to look at pecnv.sh in detail to see how these last 4 executables are called.  In the future, their interfaces will likely change as the package evolves to support BAM files generated by "bwa mem".

#Using the programs to detect CNVs (not TEs, though -- see above for our pipeline to do that).

The easiest way to use the programs will be to run the master script, _pecnv.sh_, on your data.  This script requires that the following software be in your user's path:

1. The [bwa](http://bio-bwa.sourceforge.net/) aligner.  __NOTE:__ this pipeline has only been used with bwa version 0.5.9!
2. [samtools](http://samtools.sourceforge.net/)
3. Rscript 

If you are impatient, run the script run_test_data.sh.  It will download a reference genome, two lanes of Illumina data, and then go to town.  __Ironically, this script cannot be run on a compute node of the UCI HPC.  This is because the compute nodes are networked in such a way that they cannot link back to the main node using hostname resolution, and thus the wget commands fail.__

##Running the master script:

_pecnv.sh_ takes the following options:

```
Usage: /home/krthornt/bin/pecnv.sh options
Mandatory options are:
 -i/--infile = A text file with the pairs of fastq file names listed.  One pair per line.  Must be TAB-separated
 -r/--reference = The reference file name to use for alignment.
Optional arguments are:
 -o/--outdir = output directory name for the output.  Will be created if it does not exist.
 -q/--minqual = Min. mapping quality for clustering.  Default = 30
 -m/--mismatches = Max. number alignment mismatches to allow for clustering.  Default = 3
 -g/--gaps = Max. number of alignment gaps to allow for clustering.  Default = 0
 -c/--cpu = Number of CPU to use for BWA alignment.  This should be 1/2 the number of cores you have available,
            because we align both reads from a lane at once, skipping intermediate files via process subsitution.
            Default = 32
 -s/--sortmem = Memory to be used by samtools sort.  Default = 50000000
 -a/--alnmem = Memory to be use by bwa aln step.  Default = 5000000
 -b/--bamfilebase = Prefix for bam file.  Default = pecnv_bamfile
 -u/--ulimit = MAX RAM usage for processing BAM file.  Unit is in gigabytes, e.g., 5 will be converted to 5*1024^2 bytes
Example:
/home/krthornt/bin/pecnv.sh -i readfile.txt -r reference.fa
```

###Comments on command-line options:

* The -u/--ulimit option allows the user to provide a hard RAM limit to the step where the program process_readmappings reads the BAM file.  For complex genomes with a large number of repetitively-mapping reads, the RAM usage may get quite high.  Thus, this option is provided so that the process may be killed rather than taking down the user's system.  Note that, if you use this option, you may see bizarre errors reported to stderr.  It is unlikely that these errors are actual segfaults, etc., in the program.  Rather, they are the outcome of what happens when a kill signal is sent by the system and not handled directly by the affected program.  In practice, the sorts of signals sent by ulimit violations are not always handleable, and thus the program makes no attempt to do so.
* For the -o option, . or ./ are allowed, and the output will be written to the current directory.  The -b option is used to ensure that each sample gets a unique name prefix, _e.g._  -b SAMPLEID would be a good idea, where SAMPLEID is something informative about this particular sample.

###What the script does

Starting from raw FASTQ files from a paired-end sequencing run, the steps are:

1. Align the data to a reference genome using bwa.  The alignment parameters are as described in Rogers _et al._ and Cridland _et al._.   Please note that the parameters are not the default BWA parameters.
2. The program process_readmappings reads the resulting BAM file, and collects reads in unusual mapping orientations, writing data to several output files.
3. The program bwa_mapdistance estimates the insert size distribution from the BAM file.  This is done separately from step 2 to minimize RAM use.  Also, power users can modify the work flow to separate these tasks out on a cluster.
4. Rscript is invoked to get the 99.9th quantile of the insert size distribution
5. cluster_cnv clusters the divergent, parallel, and unlinked read pairs into putative CNV calls.  The output files are described below.

###General comments on the work flow

* If the bam file exists, the script will skip the alignment step.  However, it will automatically redo the scanning and clustering steps.

##The output of the CNV clustering workflow

The primary output from _pecnv.sh_ is the following:

* $ODIR/$BAM.div.gz
* $ODIR/$BAM.par.gz
* $ODIR/$BAM.ul.gz

Where $ODIR is the value passed to the -o option and $BAM is the bam file name.

The first three files in the above list correspond to clusters of reads mapping in divergent orientation, parallel orientation, and read pairs mapping to different chromosomes, respectively.  See Figure 2 of [Cridland _et al._ 2010](http://gbe.oxfordjournals.org/content/2/83.full) for a figure representing these three mapping types.

The format of the output files is as follows:

* id = Event identification number (arb. integer)
* chrom1 = Chromosome number in reference where the first read cluster is
* coverage = Number of read pairs supporting the event
* strand1 = Strand of first read cluster.  0 = plus, 1 = minus
* start1 = Start position of first read cluster.  
* stop1 = Stop position of second read cluster.
* chrom2 = Chromosome number in reference where the second read cluster is
* start2 = Start position of second read cluster.  
* stop2 = Stop position of second read cluster. 
* reads = Pipe-separated (the pipe is the | character) list of the read pairs supporting the event.  Format is readPairName;start,stop,strand,start,stop,strand, where the last two values are for the two reads in the pair.  

FOr the above, all references to positions are with respect to a 1-offset coordinate system.  (In other words, 1 is the first position on a contig.)

For all referencs to "strand" in the above, 0 = plus, 1 = minus.

###Secondary output files

(This section is probable for power-users only.)

The program process_readmappings creates the following files:

* $ODIR/$BAM.cnv_mappings.csv.gz = reads mapping in parallel, divergent, and unlikned orientation
* $ODIR/$BAM.cnv_mappings.sam.gz = pseudo-SAM format records corresponding to the above file.  (This can be deleted, and its primary raison d'etre is debugging/double-checking.)
* $ODIR/$BAM.um_u.csv.gz = The uniquely-mapping read in a unique/repetitive read pair.
* $ODIR/$BAM.um_u.sam.gz = The pseudo-SAM corresponding to the above
* $ODIR/$BAM.um_m.csv.gz = The repetitively-mapping read in a unique/repetitive read pair.
* $ODIR/$BAM.um_m.sam.gz = The pseudo-SAM corresponding to the above.

The SAM files are really pseudo-SAM because they are quick-and dirty conversion of the binary BAM records, and have not been prettied up the way that samtools does.  However, they contain the same info in the same order.

For all of the below, strand = 0 or 1 for plus or minus, respectively.  All genomic positions start from 0, __not from 1__.

The format of $ODIR/$BAM.cnv_mappings.csv.gz is 17 columns:

1. Read name pair prefix (string). The hash symbol and remaining characters have been stripped, for reading in R.
2. Mapping quality of the read1. (int)
3. Chromosome for read1. (string)
4. Alignment start for read 1. (int)
5. Alignment stop for read 1. (int)
6. Strand for read 1. (int) 
7. Mismatches for read 1. (int)
8. Gaps in alignment of read 1 to ref. (int)
9. Read pair orientation (string).  UL = unlinked, DIV = divergent, PAR = parallel
10. Mappint quality of read 2. (int)
11. Chromosome for read 2. (string)
12. Alignment start for read 2. (int)
13. Alignment stop for read 2. (int)
14. Strand for read 2. (int) 
15. Mismatches for read 2. (int)
16. Gaps in alignment of read 2 to ref. (int)
17. Read pair orientation (string).  UL = unlinked, DIV = divergent, PAR = parallel

Yes, columns 9 and 17 are redundant, and are used to make sure that the output routines are working OK.

The format of $ODIR/$BAM.um_u.csv.gz and $ODIR/$BAM.um_m.csv.gz consist of 8 columns:

1. Read name pair prefix. (string)
2. Mapping quality of the read. (int)
3. Chromosome for read. (string)
4. Alignment start for read. (int)
5. Alignment stop for read. (int)
6. Strand for read. (int) 
7. Mismatches for read. (int)
8. Gaps in alignment of read to ref. (int)

The difference between the two is that the former contains 1 record per read name pair prefix, while the latter contains > 1.

The program bwa_mapdistance outputs 3 columns in a file called $ODIR/$BAM.mdist.gz.

1. distance (int).  This is the ISIZE field from the BAM record.
2. number (int).  The number of read pairs with that distance
3. cprob (double). The cumulative density of distance in the ECDF of ISIZES.

The program pecnv_insert_qtile is a simple Rscript that reads the output of bwa_mapdistance and outputs a value.  It takes $ODIR/$BAM.mdist.gz and a value 0 <= x <= 1 as arguments.  The distance corresponding to the x-th quantile is printed to STDOUT.  The Rscript is used internally to capture these quantile values.

##Running on pre-existing bam files

The main script, pecnv.sh expects a bamfile called PREFIX_sorted.bam, where PREFIX is the value passed to the -b/--bamfile option (see above).  If you already have an existing BAM file, but its name is not what is expected, you may make a symbolic link that has the correct name pattern, and then run the script.

#Detecting tranposable element (TE) insertions