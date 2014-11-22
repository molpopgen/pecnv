#ifndef __TECLUST_OBJECTS_HPP__
#define __TECLUST_OBJECTS_HPP__

#include <cstdint>
#include <string>
#include <utility>

struct teclust_params //Command-line parameter options
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
  std::int32_t INSERTSIZE,MDIST,MINREADS;
  /*
    Closest distance to known TE in reference.  For PHRAP output: only write if pdist || mdist > CLOSEST
  */
  int CLOSEST;
  /*
    For PHRAP output: only try to assemble novel insertions.
    Use the greedy algo of Cridland et al.?
  */
  bool novelOnly,greedy;
  teclust_params();
}; 

struct teinfo : public std::pair<std::int32_t,std::int32_t>
{
  std::int32_t start() const;
  std::int32_t stop() const;
  teinfo( std::int32_t __s, std::int32_t __st );
};

struct cluster
/*
  A cluster is a group of unique reads that suggest the presence of a TE
*/
{
  std::pair<std::int32_t,std::int32_t> positions;
  unsigned nreads;
  cluster();
  cluster(const int32_t & pos1,
	  const int32_t & pos2,
	  const unsigned & nr);
};
#endif
