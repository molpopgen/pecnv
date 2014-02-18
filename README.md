#Introduction
_pecnv_ is a collection of C++ programs used for the detection of copy-number variation (CNV) by means of paired-end mapping of short reads.  

This software contains the C++ portions of the transposable element (TE) detection pipeline descrbed in:
> Cridland, J.M., S.J. MacDonald, A.D. Long, and K.R Thornton (2013) Abundance and Distribution of Transposable Elements in Two \textit{Drosophila} QTL Mapping Resources  Molecular Biology and Evolution 30: 2311-2327, which is available [here](http://mbe.oxfordjournals.org/content/30/10/2311.full)

The full version of the TE detection pipeline, along with test data, is available from the Thornton lab [website](http://www.molpopgen.org/tepipeline/line99_example.tar.gz).  Users interested primarily in TE detection are encouraged to download that archive and study how it works.

This code has been used by the lab for detecting tandem duplications from short-read data, and the accompanying perl/shell scripts are intended to walk a user through doing such an analysis.

#Installation

##Dependencies

This software requires the following libraries to compile:

1. [libsequence](http://www.github.com/molpopgen/libsequence) - version 1.8.0 or greater
2. [boost](http://www.boost.org) - version 0.5.3 or greater

##Compilation

On most systems, simply type "make" to compile the programs.  On systems where the dependencies may be installed in funny locations, edit the Makefile accordingly.

##Installation

Either copy the binaries to somewhere in your user's path, or add the directory where you compiled them to your user's path.  To copy them into your user's ~/bin directory on a Linux machine, this will do the trick:

find . -perm -111 -type f -maxdepth 1 -exec cp "{}" ~/bin \;

#Using the programs to detect CNVs (not TEs, though -- see above for our pipeline to do that).

The easiest way to use the programs will be to run the master script, _pecnv.pl_, on your data.  This script requires that the following software be in your user's path:

1. The [bwa](http://bio-bwa.sourceforge.net/) aligner.
2. [samtools](http://samtools.sourceforge.net/)

__Note:__ the _pecnv.pl_ script uses ForkManager to attempt to make more effective use of CPU resources during several of the steps.  Make sure that this perl module is installed on your system!

#Running the master script:

_pecnv.pl_ takes the following options:

> pecnv.pl -outdir pecnv_output -minqual 30 -mismatches 3 -gaps 0 -infile (infilename) -sample 0 -cpu 32 -ref (reference_fasta_filename)

In the above line, default values are shown for each option where the exist and values in parentheses must be provided by the user.  The definition of each argument is:

1. outdir = the name of the directory to write the pipeline output.  This is ALL files, including bam files, etc.
2. minqual = minimum mapping quality to consider a read in the CNV clustering
3. mismatches = max number of mismatches to allow in a read alignment to consider it for CNV clustering
4. gaps = max number of gaps to allow in a read alignment to consider it for CNV clustering
5. infile = an input file containing names of FASTQ files.  They must be in order of left read, right read, left read, right read, for all read files from a single sample.
6. sample = an integer that will be used to identify this sample.  Must be non-negative.
7. cpu = number of CPU to use
8. ref = name of the fasta file containing the renamed reference

#The output

The (interesting) output from _pecnv.pl_ is the following:

