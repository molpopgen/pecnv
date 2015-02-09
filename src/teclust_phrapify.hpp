#ifndef __PHRAPIFY_HPP__
#define __PHRAPIFY_HPP__

#include <string>
#include <vector>
#include <teclust_objects.hpp>

void phrapify( const teclust_params & pars,
	       const std::string & clusters );

void phrapify_t( const teclust_params & pars,
		 const std::vector<std::pair<std::string,std::pair<std::uint64_t,std::uint64_t> > > & offsets,
		 const std::string & clusters );

#endif
