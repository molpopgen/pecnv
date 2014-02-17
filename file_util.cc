#include <file_util.hpp>
#include <iostream>
#include <fstream>

using namespace std;

bool file_exists(const char * fn)
/*
  not most rigorous way to do this
*/
{
  ifstream in(fn);
  if(in) 
    {
      return 1;
    }
  return 0;
}
