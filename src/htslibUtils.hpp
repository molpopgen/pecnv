#ifndef __PECNV_HTSLIB_UTILS_HPP__
#define __PECNV_HTSLIB_UTILS_HPP__

#include <string>
#include <vector>
#include <utility>
#include <cstdint>
#include <htslib/sam.h>

/*
  Return a pointer to the index of a bam file.
  nullptr is returned if no index file is found.
  The name of the index file is assumed to be bamfilename.bai
*/
hts_idx_t * get_index(const std::string & bamfilename);

/*
  Returns a set of bam file offsets pointing to the beginning and end of every chromosome to which reads are mapped.
  The return value is empty if a bam file index cannot be found.
  The name of the index file is assumed to be bamfilename.bai.
 */
std::vector<std::pair<std::string,std::pair<std::uint64_t,std::uint64_t> >> read_index(const std::string & bamfilename);

struct refrange
{
  std::uint64_t beg,end;
  std::int32_t refid;
  std::uint32_t start,stop;
  std::string refname;
  refrange( const std::uint64_t & __beg,
	    const std::uint64_t & __end,
	    const std::int32_t & __refid,
	    const std::uint32_t & __start,
	    const std::uint32_t & __stop,
	    const std::string  & __refname) : beg(__beg),end(__end),refid(__refid),
					      start(__start),stop(__stop),refname(__refname)
  {
  }
};

struct bamrange
{
  std::uint64_t beg,end;
  std::int32_t refid1,refid2;
  std::uint32_t start,stop;
  std::string refname1,refname2;
  bamrange() = default;
  bamrange( const std::vector<refrange> & ranges );
};

std::vector<bamrange> split_genome(const std::string & bamfilename,
				   const unsigned & NTHREADS);
#endif
