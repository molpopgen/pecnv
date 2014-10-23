#include <Sequence/Fasta.hpp>
#include <iostream>
#include <fstream>
#include <isbinary.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

using namespace std;
using namespace boost::iostreams;
using namespace Sequence;

int main(int argc, char ** argv)
{
  int argn=1;
  if (argc != 5)
    {
      cerr << "usage: rename_reference reference new_reference name_table seq_table\n";
      exit(0);
    }
  const char * ref=argv[argn++];
  const char * of = argv[argn++];
  const char * table = argv[argn++];
  const char * table2 = argv[argn++];

  Fasta f;
  unsigned u=0;
  ofstream tablestream(table),table2stream(table2);
  filtering_ostream newref;
  newref.push(gzip_compressor());
  newref.push(file_sink(of),ios_base::out | ios_base::binary);
  if(! isbinary(ref) )
    {
      ifstream fasta(ref);
      
      while( !fasta.eof() )
	{
	  fasta >> f >> ws;
	  string::size_type sst = f.first.find(' ');
	  string name;
	  if(sst != string::npos)
	    {
	      name=string(f.first.begin(),f.first.begin()+sst);
	    }
	  else
	    {
	      name=string(f.first.begin(),f.first.end());
	    }
	  
	  newref << '>' << u << '\n'
		 << f.second << '\n';
	  table2stream << u << '\t' << f.second << '\n';
	  tablestream <<  u++ << '\t' << name << '\n';
	}
    }
  else
    {
      filtering_istream fis;
      fis.push(gzip_decompressor());
      fis.push(file_source(ref,ios_base::in|ios_base::binary));
      while( !fis.eof() )
	{
	  fis >> f >> ws;
	  string::size_type sst = f.first.find(' ');
	  string name;
	  if(sst != string::npos)
	    {
	      name=string(f.first.begin(),f.first.begin()+sst);
	    }
	  else
	    {
	      name=string(f.first.begin(),f.first.end());
	    }
	  
	  newref << '>' << u << '\n'
		 << f.second << '\n';
	  table2stream << u << '\t' << f.second << '\n';
	  tablestream <<  u++ << '\t' << name << '\n';
	}
    }
}
