#ifndef __STRING_UNSIGNED_LOOKUP_HPP__
#define __STRING_UNSIGNED_LOOKUP_HPP__

#include <string>
#include <utility>
#include <vector>

unsigned update_lookup( std::vector<std::pair<std::string,unsigned> > * vpsu,
			unsigned * index,
			const std::string & s );

unsigned lookup_unsigned( const std::vector<std::pair<std::string,unsigned> > & vpsu,
			  const std::string & s );


std::string lookup_string( const std::vector<std::pair<std::string,unsigned> > & vpsu,
			   const unsigned & u );

#endif
