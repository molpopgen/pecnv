#ifndef __TECLUST_OBJECTS_HPP__
#define __TECLUST_OBJECTS_HPP__

#include <string>
#include <utility>

struct params //Command-line parameter options
{
  /*
    The reference TE positions, 
    output file name, 
    input bam file name, 
    file name listing the FASTQ files
    the um_u file that is the output of cluster_cnv for the sample
    the um_m files that is the output of cluster_cnv for the sample
    The output of this program run on the reference genome, if available
    Directory for writing the phrap output
  */
  std::string  reference_datafile, outfile, bamfile, readfile, umufile, ummfile, phrapdir;//,teclust_ref;
  /*
    Upper limit on insert size distribution
    Maximum distance used for matching up left and right ends of putative TE calls
    Min # of reads for phrap output, if wanted
  */
  unsigned INSERTSIZE,MDIST,MINREADS;
  params();
}; 

struct teinfo : public std::pair<unsigned,unsigned>
{
  unsigned start() const;
  unsigned stop() const;
  teinfo( unsigned __s, unsigned __st );
};

struct cluster
/*
  A cluster is a group of unique reads that suggest the presence of a TE
*/
{
  std::pair<unsigned,unsigned> positions;
  unsigned nreads;
  cluster();
  cluster(const unsigned & pos1,
	  const unsigned & pos2,
	  const unsigned & nr);
};
#endif
