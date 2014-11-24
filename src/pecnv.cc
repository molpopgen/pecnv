#include <teclust.hpp>
#include <algorithm>
#include <cstring>
#include <iostream>

using namespace std;

void usage(int status);
char ** strip_argv(int argc, char ** argv, const char * pattern);

int main( int argc, char ** argv )
{
  if( argc < 2 )
    {
      usage(0);
    }

  //What module will re run?
  if( strcmp(argv[1],"teclust") == 0 )
    {
      auto x = strip_argv(argc,argv,argv[1]);
      teclust_main(x - argv, argv);
    }
}

void usage(int status)
{

  exit(status);
}

char ** strip_argv(int argc, char ** argv, const char * pattern)
{
  return remove_if(argv,argv+argc,[&](const char * __xx) { return std::strcmp(__xx,pattern) == 0; });
}
