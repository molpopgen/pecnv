#include <Sequence/Fasta.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <file_common.hpp>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <zlib.h>

using namespace std;
using namespace boost::program_options;
using namespace boost::iostreams;
using namespace Sequence;

template<typename streamtype> 
void process( streamtype & in,
	      const char * outfilename )
{
  Fasta f;
  gzFile gzout = gzopen(outfilename,"w");
  if(gzout==NULL)
    {
      cerr << "Error: could not open "
	   << outfilename 
	   << " for writing on line "
	   << __LINE__ << " of " << __FILE__ << '\n';
      exit(1);
    }
  ostringstream buffer;
  while(!in.eof())
    {
      in >> f >> ws;
      auto __w = find_if(f.first.begin(),f.first.end(),[](const char & ch) { return isspace(ch); });
      buffer << string(f.first.begin(),__w) << '\t' << f.second.size();
      if(!in.eof()) buffer << '\t';
    }
  buffer << '\n';
  if ( gzwrite(gzout,buffer.str().c_str(),buffer.str().size()) <= 0 )
    {
      cerr << "Error: gzwrite error on line " << __LINE__
	   << " of " << __FILE__ << '\n';
      exit(1);
    }
  gzclose(gzout);
}

int mkgenome_main( int argc, char ** argv )
{
  string infile(argv[1]),outfile(argv[2]);
  if( file_exists(infile.c_str()) )
    {
      if(is_gzip(infile.c_str()))
	{
	  cerr << "gz\n";
	  filtering_istream in;
	  in.push(gzip_decompressor());
	  in.push(file_source(infile.c_str()),ios_base::in|ios_base::binary);
	  process(in,outfile.c_str());
	}
      else
	{
	  cerr << "ascii " << argv[1] << '\n';
	  ifstream in(infile.c_str());
	  process(in,outfile.c_str());
	}
    }
  return 0;
}


