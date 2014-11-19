#ifndef __PECNV_INTERMEDIATEIO_HPP__
#define __PECNV_INTERMEDIATEIO_HPP__

#include <Sequence/bamrecord.hpp>
#include <zlib.h>
#include <cstdint>
#include <string>
#include <utility>

int gzwriteCstr( gzFile, const std::string & );
std::pair<std::string,int> gzreadCstr( gzFile in);

struct alnInfo
{
  std::int32_t start,stop;
  std::int8_t mapq,strand;
  std::int16_t mm,ngap; //no mismatchs, gaps, resp.

  alnInfo( const Sequence::bamrecord & b ); //construct from an alignment
  alnInfo( gzFile in );//construct via input from a gzFile;
  alnInfo( const int32_t &,
	   const int32_t &,
	   const int32_t &,
	   const int32_t &,
	   const uint32_t &,
	   const uint32_t & ); //construct from raw numbers
  int write( gzFile in);
};

#endif

