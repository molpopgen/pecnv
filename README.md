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

find . -perm -111 -type f -maxdepth 1 -exec cp "{}" ~/bin \

#Using the programs to detect CNVs (not TEs, though -- see above for our pipeline to do that).

The easiest way to use the programs will be to run the master script, _pecnv.pl_, on your data.  This script requires that the following software be in your user's path:

1. The [bwa](http://bio-bwa.sourceforge.net/) aligner.
2. [samtools](http://samtools.sourceforge.net/)