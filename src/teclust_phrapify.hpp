#ifndef __PHRAPIFY_HPP__
#define __PHRAPIFY_HPP__

#include <string>
#include <vector>
#include <teclust_objects.hpp>
#include <htslibUtils.hpp>

void phrapify( const teclust_params & pars,
	       const std::string & clusters );

//1 thread per sequence in the reference
void phrapify_t( const teclust_params & pars,
		 const std::vector<std::pair<std::string,std::pair<std::uint64_t,std::uint64_t> > > & offsets,
		 const std::string & clusters );

//The references is divided up into pars.NTHREADS equal-sized chunks
void phrapify_t_v2( const teclust_params & pars,
		    const std::vector<bamrange> & branges,
		    const std::string & clusters );

#endif
