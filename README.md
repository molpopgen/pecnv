#Introduction
_pecnv_ is a collection of C++ programs used for the detection of copy-number variation (CNV) by means of paired-end mapping of short reads.  

This pipeline corresponds to the clustering algorithm described and used in

*  __Rogers, R. L.__, __J. M. Cridland__, __L. Shao__, T. T. Hu, P. Andolfatto and K. R. Thornton (2014) Landscape of standing variation for tandem duplications in _Drosophila yakuba_ and _Drosophila simulans_.  Molecular Biology and Evolution __31__: 1750-1766 PMID 24710518, [Manuscript](http://mbe.oxfordjournals.org/content/31/7/1750.abstract.html)

A preprint of the manuscript is available on [arXiv](http://arxiv.org/abs/1401.7371).

This software also contains the C++ portions of the transposable element (TE) detection pipeline descrbed in:

* __Cridland, J.M.__, S.J. MacDonald, A.D. Long, and K.R Thornton (2013) Abundance and Distribution of Transposable Elements in Two _Drosophila_ QTL Mapping Resources  Molecular Biology and Evolution __30__: 2311-2327. PMID 23883524 [Manuscript](http://mbe.oxfordjournals.org/content/30/10/2311.full)

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

##Dependencies

This software requires the following libraries to compile:

1. [libsequence](http://www.github.com/molpopgen/libsequence) - version 1.8.3 or greater
2. [zlib](http://zlib.net) (This is a also a libsequence dependency)
3. [htslib](http://htslib.org)

__Please make sure that your installation of libsequence makes use of htslib.  That requires that htslib be installed prior to installing libsequence, so that the latter will compile features depending on the former.__

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

##Using the programs to detect CNVs (not TEs, though -- see above for our pipeline to do that).

The easiest way to use the programs will be to run the master script, _pecnv.sh_, on your data.  This script requires that the following software be in your user's path:

1. The [bwa](http://bio-bwa.sourceforge.net/) aligner.  __NOTE:__ this pipeline has only been used with bwa version 0.5.9!
2. [samtools](http://samtools.sourceforge.net/)
3. Rscript 

If you are impatient, run the script run _ test _ data.sh.  It will download a reference genome, two lanes of Illumina data, and then go to town.  __Ironically, this script cannot be run on a compute node of the UCI HPC.  This is because the compute nodes are networked in such a way that they cannot link back to the main node using hostname resolution, and thus the wget commands fail.__

###Running the master script:

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

####Comments on command-line options:

* The -u/--ulimit option allows the user to provide a hard RAM limit to the step where the program process_readmappings reads the BAM file.  For complex genomes with a large number of repetitively-mapping reads, the RAM usage may get quite high.  Thus, this option is provided so that the process may be killed rather than taking down the user's system.  Note that, if you use this option, you may see bizarre errors reported to stderr.  It is unlikely that these errors are actual segfaults, etc., in the program.  Rather, they are the outcome of what happens when a kill signal is sent by the system and not handled directly by the affected program.  In practice, the sorts of signals sent by ulimit violations are not always handleable, and thus the program makes no attempt to do so.
* For the -o option, . or ./ are allowed, and the output will be written to the current directory.  The -b option is used to ensure that each sample gets a unique name prefix, _e.g._  -b SAMPLEID would be a good idea, where SAMPLEID is something informative about this particular sample.

####What the script does

Starting from raw FASTQ files from a paired-end sequencing run, the steps are:

1. Align the data to a reference genome using bwa.  The alignment parameters are as described in Rogers _et al._ and Cridland _et al._.   Please note that the parameters are not the default BWA parameters.
2. The program process_readmappings reads the resulting BAM file, and collects reads in unusual mapping orientations, writing data to several output files.
3. The program bwa_mapdistance estimates the insert size distribution from the BAM file.  This is done separately from step 2 to minimize RAM use.  Also, power users can modify the work flow to separate these tasks out on a cluster.
4. Rscript is invoked to get the 99.9th quantile of the insert size distribution
5. cluster_cnv clusters the divergent, parallel, and unlinked read pairs into putative CNV calls.  The output files are described below.

####General comments on the work flow

* If the bam file exists, the script will skip the alignment step.  However, it will automatically redo the scanning and clustering steps.

###The output

The (interesting) output from _pecnv.sh_ is the following:

* $ODIR/$BAM.div.gz
* $ODIR/$BAM.par.gz
* $ODIR/$BAM.ul.gz

Where $ODIR is the value passed to the -o option and $BAM the value passed to the -b option.  

The above correspond to clusters of reads mapping in divergent orientation, parallel orientation, and read pairs mapping to different chromosomes, respectively.  See Figure 2 of [Cridland _et al._ 2010](http://gbe.oxfordjournals.org/content/2/83.full) for a figure representing these three mapping types.

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

###Running on pre-existing bam files

The main script, pecnv.sh expects a bamfile called PREFIX_sorted.bam, where PREFIX is the value passed to the -b/--bamfile option (see above).  If you already have an existing BAM file, but its name is not what is expected, you may make a symbolic link that has the correct name pattern, and then run the script.