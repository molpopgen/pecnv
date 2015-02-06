#ifndef __TECLUST_SCAN_BAMFILE_HPP__
#define __TECLUST_SCAN_BAMFILE_HPP__

#include <teclust_objects.hpp>
#include <unordered_set>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <cstdint>
#include <htslib/sam.h>
using refTEcont = std::map<std::string,std::vector<teinfo> >;

void scan_bamfile(const teclust_params & p,
		  const refTEcont & refTEs,
		  std::unordered_set<std::string> * readPairs,
		  std::map<std::string,std::vector< std::pair<std::int32_t,std::int8_t> > > * data,
		  hts_idx_t * idx);

#endif
