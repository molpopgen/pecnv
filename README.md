#Introduction
_pecnv_ is a collection of C++ programs used for the detection of copy-number variation (CNV) by means of paired-end mapping of short reads.  

This pipeline corresponds to the clustering algorithm described and used in

*  Rogers, R. L.,J. M. Cridland, L. Shao, T. T. Hu, P. Andolfatto and K. R. Thornton (2014) Landscape of standing variation for tandem duplications in _Drosophila yakuba_ and _Drosophila simulans_.  Molecular Biology and Evolution __31__: 1750-1766 PMID 24710518, [Manuscript](http://mbe.oxfordjournals.org/content/31/7/1750.abstract.html)

This package also contains a C++ implementation of the transposable element (TE) detection pipeline descrbed in:

* Cridland, J.M., S.J. MacDonald, A.D. Long, and K.R Thornton (2013) Abundance and Distribution of Transposable Elements in Two _Drosophila_ QTL Mapping Resources  Molecular Biology and Evolution __30__: 2311-2327. PMID 23883524 [Manuscript](http://mbe.oxfordjournals.org/content/30/10/2311.full)

#A note on version numbers

1.  Release 0.1.0 = the precise version used in Rogers et al. and Cridland et al.
2.  Release 0.1.1 = a bug fix resolving issue number 1 regarding inconsistent SAM flags within read pairs.  This bug does not affect the conclusions of the previous work because it affects very few reads that make it into the final clustering steps.
3.  Release 0.1.2 = No more need to rename references or read names.  This shaves a load of time off of the workflow.
4.  Release 0.1.3 = re-implementation of process_readmappings to read directly from BAM files.  No more need to make a BAM file sorted by read name.  This shaves hours off of the pipeline.
5.  Release 0.1.4 = integration of the workflow for calling transposable-element insertions from PE mapping data using the Cridland _et al._ approach.  Intermediate files are now written in a binary format for more efficient downstream-processing.  This version depends on the boost libraries for parsing command-line arguments

##Tentative plan for future releases

1. Bugfixes, better documentation, etc.
2. Support for additional output files in BED format, for compatibility with other tools.

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

These run-time dependencies must all be in a user's path.  Further, the system must be set up like a "normal" *nix machine so that a script may call /usr/bin/env XXX to find the path to XXX on the [shebang](http://en.wikipedia.org/wiki/Shebang_(Unix)) line of a script.

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

* pecnv.sh is an executable bash script that runs the CNV calling pipeline on paired-end short-read Illumina data (see next section)
*teclust is a C++ program implementing the TE detection method of Cridland _et al._ (2013) PMID 23883524.  This program is explained in more detail below.

The following executables are installed and called by pecnv.sh:
* process_readmappings is a C++ program that collects read pairs in unxepected orientations (divergent, parallel, unlinked, and unique/multi) from a BAM file
* bwa_mapdistance is a C++ program that summarizes the ECDF of the insert size distribution based on uniquely-mapped paired reads in a BAM file
* pecnv_insert_qtile is an executable Rscript that processes the output of bwa_mapdistance to estimate the upper quantile of the insert size distribtution
* cluster_cnv is a C++ program that clusters the output of process_readpmappings into putative CNV calls

"Power users" are encouraged to look at pecnv.sh in detail to see how these last 4 executables are called.  In the future, their interfaces will likely change as the package evolves to support BAM files generated by "bwa mem".

#Using the programs to detect CNVs

The easiest way to use the programs will be to run the master script, _pecnv.sh_, on your data.  This script requires that the following software be in your user's path:

1. The [bwa](http://bio-bwa.sourceforge.net/) aligner.  __NOTE:__ this pipeline has only been used with bwa version 0.5.9!
2. [samtools](http://samtools.sourceforge.net/)
3. Rscript 

If you are impatient, run the script run_test_data.sh, found in the test subdirectory of the source code repository.  It will download a reference genome, two lanes of Illumina data, and then go to town.  __Ironically, this script cannot be run on a compute node of the UCI HPC.  This is because the compute nodes are networked in such a way that they cannot link back to the main node using hostname resolution, and thus the wget commands fail.__

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
* $ODIR/$BAM.unl.gz

Where $ODIR is the value passed to the -o option and $BAM is the bam file name.

The first three files in the above list correspond to clusters of reads mapping in divergent orientation, parallel orientation, and read pairs mapping to different chromosomes ("unlinked"), respectively.  See Figure 2 of [Cridland _et al._ 2010](http://gbe.oxfordjournals.org/content/2/83.full) for a figure representing these three mapping types.

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

###What does the output mean?

If you are working in a system with an incomplete genome, then reads mapping to differnt contigs should treated with some caution, as you cannot assume that you know the true mapping relationship of those reads.

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

The format of $ODIR/$BAM.cnv_mappings.csv.gz a gzipped file binary-format records describing unusual read mappings.  Note that this file is not human-readable.  The format of a single record is:

1. Read name pair prefix (\0-terminated C string). The hash symbol and remaining characters have been stripped.
2. Chromosome name for read1 (\0-terminated C string).
3. Chromosome name for read2 (\0-terminated C string).
4. Read pair orientation (char[3]). UNL = unlinked, DIV = divergent, PAR = parallel
5. Alignment start position for read 1 (int32_t).
6. Alignment stop position for read 1 (int32_t).
7. Mapping quality for read 1 (int8_t)
8. Strand for read 1 (int8_t).
9. Mismatches for read 1 (int16_t)
10. Alignment gaps for read 1 (int16_t)
11. Alignment start position for read 2 (int32_t).
12. Alignment stop position for read 2 (int32_t).
13. Mapping quality for read 2 (int8_t)
14. Strand for read 2 (int8_t).
15. Mismatches for read 2 (int16_t)
16. Alignment gaps for read 2 (int16_t)

The values in parentheses correspond to the data types using in C/C++.  intX_t refers to a signed integer guaranteed to be exactly X bits in size.  Your system defines these types in <stdint.h> (C) and <cstdint> (C++11).  This project manages the IO for these files using the routines defined in the file __intermediateIO.hpp__. A \0-terminated C string is a "char *" including the C-language string termination character \0. These data fields are read by repeated calls to gzgetc until the \0 is found.

The format of $ODIR/$BAM.um_u.csv.gz and $ODIR/$BAM.um_m.csv.gz are again binary-encoded records in the following format:

1. Read name pair prefix (\0-terminated C string). The hash symbol and remaining characters have been stripped.
2. Chromosome for read. (\0-terminated C string)
3. Alignment start for read. (int32_t)
4. Alignment stop for read. (int32_t)
5. Mapping quality of the read. (int8_t)
6. Strand for read. (int8_t) 
7. Mismatches for read. (int16_t)
8. Gaps in alignment of read to ref. (int16_t)

The difference between the two is that the former contains 1 record per read name pair prefix, while the latter contains > 1.

The program bwa_mapdistance outputs 3 columns in a file called $ODIR/$BAM.mdist.gz.

1. distance (int).  This is the ISIZE field from the BAM record.
2. number (int).  The number of read pairs with that distance
3. cprob (double). The cumulative density of distance in the ECDF of ISIZES.

The program pecnv_insert_qtile is a simple Rscript that reads the output of bwa_mapdistance and outputs a value.  It takes $ODIR/$BAM.mdist.gz and a value 0 <= x <= 1 as arguments.  The distance corresponding to the x-th quantile is printed to STDOUT.  The Rscript is used internally to capture these quantile values.

##Running on pre-existing bam files

The main script, pecnv.sh expects a bamfile called PREFIX_sorted.bam, where PREFIX is the value passed to the -b/--bamfile option (see above).  If you already have an existing BAM file, but its name is not what is expected, you may make a symbolic link that has the correct name pattern, and then run the script.

#Detecting tranposable element (TE) insertions

The Cridland et al. method is implemented in the program teclust, with some additional features.  The program is able to make use of TE annotations from an existing reference genome, if they are available.

The minimal input for this program is the output of the pecnv.sh pipeline, which needs to be run first.  Specifically, the "um_u" and "um_m" files described above are needed.

There is a test script in the test subdirectory of the source code repository.  See the README.md in that directory for more details.

The usage information for the program is:

```
Cluster reads into putative transposable element calls.
Usage: tclust -h to see help:
  -h [ --help ]                Produce help message
  -b [ --bamfile ] arg         BAM file name (optional)
  -t [ --tepos ] arg           File containing positions of TEs in reference 
                               genome (optional)
  -o [ --outfile ] arg         Output file name for clusters (required)
  -u [ --umu ] arg             The um_u output file for the sample generated by
                               process_readmappings (required)
  -m [ --umm ] arg             The um_u output file for the sample generated by
                               process_readmappings (required)
  -i [ --isize ] arg           Upper limit on insert size distribution, e.g. 
                               from bwa_mapdistance. (required)
  -M [ --mdist ] arg (=1000)   Max. distance for joining up left and right ends
                               of putative TEs (required)
  -p [ --phrapdir ] arg        Name of a directory to put input files for de 
                               novo assembly of putatitve TE insertions using 
                               phrap. If the directory does not exist, it will 
                               be created. (optional)
  -r [ --minreads ] arg (=3)   Min. number of reads in a cluster for writing 
                               input files for phrap. (optional)
  -c [ --closestTE ] arg (=-1) For phrap output, only consider events >= c bp 
                               away from closest TE in the reference. 
                               (optional)
  --ummHitTE                   When processing the um_u/um_m files, only 
                               consider reads where the M read hits a known TE.
                                 This makes --tepos/-t a required option. 
                               (optional)
  -a [ --allEvents ]           For phrap output: write files for all events. 
                               Default is only to write files for putative 
                               novel insertions
```

##Various ways to do it

###No annotated TEs in your reference.

The simplest command line possible is:

```
#Use the 99th quantile of mapping distances
teclust -u umuFile -m ummFile -o outfile.gz -i `pecnv_insert_qtile mdistfile 0.99`
```

The output of this analysis will be a set of clusters defined by a run of uniquely-aligning reads whose partners to not map uniquely.

###Annotated TEs in the reference genome

If you work with a system where you can provide a list of TE positions in the reference genome, you can do things differently.

First, you can require that the M reads in the umm file overlap with an annotated TE position:

```
teclust -u umuFile -m ummFile -o outfile.gz -i `pecnv_insert_qtile mdistfile 0.99` -t tefile --ummHitTE
```

The output is now restricted to clusters defined the the overlap of the non-unique read with an annotated TE.

The format of "tefile" is 3 columns:

* chromosome name (string).  These names MUST be the same as what is found in your BAM files.
* start (int)
* stop (int)

Start and stop have ranges (1,chromlen).  Start must always be less than stop, and the prorgram __does not__ check this!  Further, the program must end with a newline character, else the last line will not be read in.

The "tefile" may be either plain text or gzipped text.

For the record, the following commands are equivalent:

```
#Pass teclust a list of TE positions, but don't require that the contents of the umm file overlap those positions
teclust -u umuFile -m ummFile -o outfile.gz -i `pecnv_insert_qtile mdistfile 0.99` -t tefile
#Run the program as you would for an unannotated genome
teclust -u umuFile -m ummFile -o outfile.gz -i `pecnv_insert_qtile mdistfile 0.99` 
```

The reason is that you're not asking the program to use the info in "tefile".

####Comments on the TE input file

For our published work, the TE input file was based on identifying regions of the genome with significant homology matches to at least 75% of the canonincal sequence for a known TE, based on the file transposon_sequence_set.embl.txt.gz from FlyBase.  That approach worked rather well.

In the test data, I show how to quickly make a TE input file using a quick combo of blast and bedtools.  I'm able to find 99 out of the 100 simulated insertions in the test data using this file (see TEvariationVariationsOnTheme.sh in the test directory of this repo for command lines using my new input file with test data), compared to 85 out of 100 using the published pipeline.

###The Cridland _et al.__ 2013 approach

To implement the procedure from this paper:

```
teclust -u umuFile -m ummFile -o outfile.gz -i `pecnv_insert_qtile mdistfile 0.99` -t tefile -b bamfile
```

The program will scan through the bamfile and look at all reads mapping to a known TE (unique or not), and use that info to augment the info in the umm/umu files.  This aids the detection of TE loci shared betweeen the sample and the reference because there are likely to be reads mapping uniquely to that specific element.

You may add the --ummHitTE option to change the protocol to require that umm data hit annotated TEs.  However, that is not what we did in the papers.

###Extracting reads for _de novo_ assembly

If you specify a directory name with the --prhapdir option, teclust will make fasta and fasta.qual files that you may run through phrap.

```
teclust -u umuFile -m ummFile -o outfile.gz -i `pecnv_insert_qtile mdistfile 0.99` -t tefile -b bamfile --phrapdir phrapdir
```

If phrapdir does not exist, the program will create it for you.  Please be careful here: if phrapdir does exist (say from an earlier run of teclust), its contents will not be affected unless the program tries to write a new file with the same name as an existing file.  In that case, the existing file will be over-written.

The file names in phrapdir have the format chrom.start.stop.[left|right].fasta and chrom.start.stop.[left|right].fasta.qual.  The left or right corresponds to the left or right cluster of a putative event.

The behavior of collecting reads for _de novo_ assembly are affected by the following teclust options

* -m/--minreads Only collect reads for events with at least 'm' reads for both the left AND the right cluster (this is based on the "pin" and "min" values -- see below in the section documenting the output file format)
* -c/--closestTE Only collect reads for evetns at least 'c' base pairs away from an annotated TE.  The default value of -1 means all events will pass this test.  Using -c 1000 will require than an event be at least 1kb away from an annotated TE.  This is enforced by a comparison to the mdist and pdist fields in the output file (see below).  If there is no annotated TE information, using -c 0 or greater will result in no events being selected for assembly, because the mdist and pdist values will be set to -1 in this case (again, see below)
* -a/--allEvents Pull reads for all events (_e.g._, novel events and TE insertions shared with the reference).  The default is putative novel insertions.  This is based on checking the "pin" and "min" fields for non-zero values in the output file.

####Automating the _de novo_ assembly.

If you want to assemble the events using [phrap](http://www.phrap.org), you need a version of phrap >= 1.090518, which is what we used in the various papers.  When we were doing the work for those papers, this version of phrap had to be specially requested from the authors and the "regular" distributed version would crash when run on short-read data.  Please be aware of what version you end up with, and be sure to get the right one.

For many users, the easiest way to automate running phrap will be via [GNU Parallel](http://www.gnu.org/software/parallel/), which is a command-line tool for automating batch jobs.  I'd suggest that you look at the [tutorial](http://www.gnu.org/software/parallel/parallel_tutorial.html) for parallel in detail, as it is quite helpful.

Here's one way to do it:

```
find phrapdir -name "*.fasta" | parallel phrap --jobs 8 --timeout 300 "phrap {} -vector_bound 0 -forcelevel 10 -minscore 10 -minmatch 10 -new_ace 2> {}.stderr > {}.stdout"
```

The command shown above does the following:

1. Use [GNU find](http://linuxcommand.org/man_pages/find1.html) to get a list of all the fasta files make by teclust
2. Use parallel to run up to 8 assemblies at a time (adjust this number for your own system's resources) and kill any job taking more than 300 seconds.
3. Run phrap with the parameters -vector_bound 0 -forcelevel 10 -minscore 10 -minmatch 10 -new_ace on each of the fasta files

The {} evaluate to each fasta file found in phrapdir.  Thus, the stderr and stdout streams from phrap are stored separately for each file.

Why the timeout?  For real data, the fasta files can get massive and the assemblies can take a very long time.  In the future, we will likely add an option to limit the number of reads in a file.  But for now, we have learned through experience that if an assembly doesn't go rather quickly, it takes a very long time.  You may have to tune the timeout duration for your system/data/etc.

You can then blast the contigs against a database of known TE sequences and use Julie Cridland's [TE annotation pipeline](https://github.com/ThorntonLab/Cridland2013AnnotPipeline) to process the output.  (The annotation methods may get merged into this project in the future.)

"Power users" may prefer array jobs on a cluster to automate the assemblies, blasting, etc.

##The output file

This section documents the contents of the file created by the --outfile/-o option to teclust.

###Background

Fundamentally, teclust groups together reads on the same strand that fall into the following categories:

* uniquely-mapping and their partners map to multiple genomic locations in the reference
* uniquely-mapping and their partners map to multiple genomic locations in the reference, at least one of which hits an annotated TE
* map near, but not in, an annotated TE, and their partners map within an annotated TE

Thus, a cluster of these mapped reads on the plus strand are assumed to be 5' of a putative TE insertion that we would guess is close by.  Likewise, a cluster of reads on the minus strand are assumed to be 3' of a putative insertion.  Further, if we obtain a large number of reads in plus- and minus- strand clusters, we infer that we're flanking the insertion site, and the distance between the ends of these clusters is an estimate of the insertion site.

Keeping this picture in mind will hopefully make the next section easier to understand.

###The file contents

The output file is gzipped text and contains a header line for easy processing using tools like R.  The columns are:

1. chromo (string) = the chromosome
2. nplus (int) = the number of read pairs on the left-hand side/plus strand of the event.
3. minus (int) = the number of read pairs on the right-hand side/minus strand of the event.  In total, columns 2 and 3 are the "coverage" in support of this event
4. pfirst (int) = the first position of the left-hand cluster
5. plast (int) = the last position of the left-hand cluster
6. pdist (int) = the distance from the left-hand cluster to the first annotated TE 3' of the cluster
7. pin (bool) = 0 of the interval pfirst to plast does __not__ overlap an annotated TE, 1 if it __does__
8. mfirst (int) = the first position of the right-hand cluster
9. mlast (int) = the last position of the right-hand cluster
10. mdist (int) = the distance from the right-hand cluster to the first annotated TE 5' of the cluster
11. min (bool) = 0 of the interval pfirst to plast does __not__ overlap an annotated TE, 1 if it __does__

The following comments apply to this file:

* positions start from 1
* The pfirst/plast and mfirst/mlast intervals are based on the alignment start positions of reads.  In other words, they correpond to the run of positions given by what would be in the fourth column of a SAM record.
* A value of -1 is used to represent "not applicable".  Old versions of this code used to output NA, which was fine when the file was read using R, but problematic when read using C/C++, etc.
* Values of -1 may occur if: nplus and/or nminus equal 0.  The --tefile/-t option is not used, therfore pdist,pin,mdist, and min cannot be known.

