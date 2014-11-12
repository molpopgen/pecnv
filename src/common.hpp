#ifndef __PECNV_COMMON_HPP__
#define __PECNV_COMMON_HPP__

#include <string>
#include <vector>
#include <utility>

//Turn a read name from XXX#something to XXX
std::string editRname( const std::string & readname );
//unsigned alen(const std::vector< std::pair<char,unsigned> > & cigar_data);
//unsigned mm(const unsigned & nm,  const std::vector< std::pair<char, unsigned> > & cigar_data);
//unsigned ngaps(const std::vector< std::pair<char,unsigned> > & cigar_data);
//std::vector<std::pair<char,unsigned> > parse_cigar(const std::string & cigar);
#endif 
