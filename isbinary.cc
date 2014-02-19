#include <isbinary.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
//#include <iostream>
//#include <fstream>
//#include <cctype>

using namespace std;
using namespace boost::iostreams;

bool isbz2( const char * filename )
{
  filtering_istream gzin;
  gzin.push(bzip2_decompressor());
  gzin.push(file_source(filename),ios_base::in|ios_base::binary);

  char c;
  try
    {
      gzin >> c;
    }
  catch ( bzip2_error & e )
    {
      gzin.pop();
      return false;
    }
  return true;
}

bool isbinary( const char * filename )
{
  filtering_istream gzin;
  gzin.push(gzip_decompressor());
  gzin.push(file_source(filename),ios_base::in|ios_base::binary);

  char c;
  try
    {
      gzin >> c;
    }
  catch ( gzip_error & e )
    {
      std::cerr << "caught exception!\n";
      gzin.pop();
      return false;
    }
  return true;
}
