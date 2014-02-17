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

  //    ifstream in(filename);
//     while( (c = in.get()) )
//       {
// 	if (! isprint(c))
// 	  {
// 	    return true;
// 	  }
// 	else
// 	  {
// 	    return false;
// 	  }
//       }

bool isbz2( const char * filename )
{
  filtering_istream gzin;
  gzin.push(bzip2_decompressor());
  gzin.push(file_source(filename));


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
  gzin.push(file_source(filename));


  char c;
  try
    {
      gzin >> c;
    }
  catch ( gzip_error & e )
    {
      gzin.pop();
      return false;
    }
  return true;
}
