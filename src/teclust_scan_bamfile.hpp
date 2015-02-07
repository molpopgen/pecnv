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
		  const hts_idx_t * idx);

void scan_bamfile_t(const unsigned & task_idx,
		    const int64_t & begin,
		    const int64_t & end,
		    const teclust_params & p,
		    const refTEcont & refTEs,
		    std::unordered_set<std::string> * readPairs,
		    std::vector< std::map<std::string,std::vector<std::pair<std::int32_t,std::int8_t> > > > & data_v,
		    const hts_idx_t * idx);
#endif
