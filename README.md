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

##Tentative plan for future releases

1.  0.1.3 will re-implement process_readmappings to read directly from BAM files.  Further, the nead for a readname-sorted BAM file will go away.  These two improvements will further speed up the pipeline.
2.  0.1.4 will integrate the workflow for calling transposable-element insertions from PE mapping data using the Cridland _et al._ approach

#Installation

##Dependencies

This software requires the following libraries to compile:

1. [libsequence](http://www.github.com/molpopgen/libsequence) - version 1.8.3 or greater
2. [zlib](http://zlib.net) (This is a also a libsequence dependency)

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
Example:
/home/krthornt/bin/pecnv.sh -i readfile.txt -r reference.fa
```

For the -o option, . or ./ are allowed, and the output will be written to the current directory.  The -b option is used to ensure that each sample gets a unique name prefix, _e.g._  -b SAMPLEID would be a good idea, where SAMPLEID is something informative about this particular sample.

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

```
#Process the (readname-sorted) bam file, and collect abnormal PE mappings:
samtools view -f 1 readsorted_bam.bam | process_readmappings cnv_mappings um mdist.gz
#Get the 99.95th quantile of fragment sizes based on unique/unique read pairs:
Rscript -e "x=read.table(\"mdist.gz\",header=T);z=which(x\$cprob >= 0.999);y=x\$distance[z[1]];write(y,\"mquant.txt\")"
#Store that value in a shell variable:
MAXDIST=`head -n 1 mquant.txt`
#Cluster the reads:
cluster_cnv minqual max_mismatches max_gaps $MAXDIST div_clusters.gz par_clusters.gz ul_clusters.gz cnv_mappings.csv.gz
```