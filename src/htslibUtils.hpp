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

#endif
