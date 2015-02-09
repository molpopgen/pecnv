#ifndef __PECNV_HTSLIB_UTILS_HPP__
#define __PECNV_HTSLIB_UTILS_HPP__

#include <string>
#include <vector>
#include <utility>
#include <cstdint>
#include <htslib/sam.h>

hts_idx_t * get_index(const std::string & bamfilename);
std::vector<std::pair<std::string,std::pair<std::uint64_t,std::uint64_t> >> read_index(const std::string & bamfilename);

#endif
