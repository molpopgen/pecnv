#ifndef __TECLUST_SCAN_BAMFILE_HPP__
#define __TECLUST_SCAN_BAMFILE_HPP__

#include <teclust_objects.hpp>
#include <unordered_set>
#include <map>
#include <string>
#include <utility>
#include <vector>

using refTEcont = std::map<std::string,std::vector<teinfo> >;

std::unordered_set<std::string> scan_bamfile(const params & p,
					     const refTEcont & refTEs,
					     std::map<std::string,std::vector< std::pair<unsigned,unsigned> > > * data);

#endif
