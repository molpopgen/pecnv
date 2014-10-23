#ifndef __BWA_UTIL_HPP__
#define __BWA_UTIL_HPP__

#include <string>
#include <Sequence/samrecord.hpp>
/*
  void getLanePair(unsigned * lane,
  unsigned * rpair,
  const std::string & name);
*/

bool hasXT(const Sequence::samrecord & r,const std::string & value = "U");

#endif


