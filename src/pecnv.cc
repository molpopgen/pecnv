#include <pecnv_version.hpp>
#include <teclust.hpp>
#include <process_readmappings.hpp>
#include <cluster_cnv.hpp>
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
  if ( strcmp(argv[1],"help") == 0 )
    {
      cout << "pecnv version " << PECNV_VERSION << '\n';
      usage(0);
    }
  if ( strcmp(argv[1],"version") == 0 )
    {
      usage(0);
    }
  else if ( strcmp(argv[1],"process") == 0 )
    {
      auto x = strip_argv(argc,argv,argv[1]);
      process_readmappings_main(x - argv, argv);
    }
  else if( strcmp(argv[1],"cnvclust") == 0 )
    {
      auto x = strip_argv(argc,argv,argv[1]);
      cluster_cnv_main(x - argv, argv);
    }
  else if( strcmp(argv[1],"teclust") == 0 )
    {
      auto x = strip_argv(argc,argv,argv[1]);
      teclust_main(x - argv, argv);
    }
}

void usage(int status)
{
  cerr << "pecnv version " << PECNV_VERSION << '\n';
  exit(status);
}

char ** strip_argv(int argc, char ** argv, const char * pattern)
{
  return remove_if(argv,argv+argc,[&](const char * __xx) { return std::strcmp(__xx,pattern) == 0; });
}
