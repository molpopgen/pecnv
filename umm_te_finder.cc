/*
  Copyright 2010 Kevin Thornton, UCI

  This program generates the input for teclust
*/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <string>
#include <functional>
#include <algorithm>
#include <isbinary.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

using namespace std;
using namespace boost::iostreams;

typedef map<string,unsigned> maptype;
struct teinfo
{
  unsigned chrom,start,stop;
  teinfo( const unsigned & c, const unsigned & st, const unsigned & sto) :
  chrom(c),start(st),stop(sto)
  {
  }
};

//the find criterion
struct iste : public binary_function<teinfo,teinfo,bool>
{
  inline bool operator()(const teinfo & lhs, const teinfo & rhs) const
  {
    bool ischrom = ( lhs.chrom == rhs.chrom ) ? true : false;
    if( !ischrom ) return false;
    bool overlap = ( (rhs.start >= lhs.start && rhs.start <= lhs.stop) || 
		     (rhs.stop >= lhs.start && rhs.stop <= lhs.stop) );
    return overlap;
  }
};


template<typename streamtype>
void read_data( maptype & umm,
		vector<teinfo> & refdata,
		streamtype & in,
		const unsigned & minmqual)
{
  //unsigned line,lane,pair,read,
  string readname;
  unsigned mq;
  unsigned chrom,start,stop;
  string temp;
  while( ! in.eof() )
    {
      //in >> line >> lane >> pair >> read >> mq
      in >> readname >> mq
	 >> chrom >> start >> stop;
      getline(in,temp);
      in >> ws;

      if( mq >= minmqual ) //filter on mapping quality
	{
	  //is this guy a TE?
	  if(!refdata.empty())
	    {
	      vector<teinfo>::iterator itr = find_if(refdata.begin(),refdata.end(),
						     bind2nd(iste(),
							     teinfo(chrom,start,stop)));
	      if ( itr != refdata.end() )
		{
		  //umm.insert(pair);
		  umm[readname]=1;
		}
	    }
	  else
	    {
	      //umm.insert(pair);
	      umm[readname]=1;
	    }
	}
    }
}

template<typename istreamtype,
	 typename ostreamtype>
void read_umu( istreamtype & in,
	       ostreamtype & out,
	       maptype & pairs)
{
  //unsigned line,lane,pair,read,mq;
  string readname;
  unsigned mq;
  unsigned chrom,start,stop,strand;
  string temp;
  unsigned nfound = 0;
  while( ! in.eof() )
    {
      //in >> line >> lane >> pair >> read >> mq
      in >> readname >> mq
	 >> chrom >> start >> stop >> strand;
      getline(in,temp);
      in >> ws;
      //if( find(pairs.begin(),pairs.end(),pair) != pairs.end() )
      //if(pairs.find(pair) != pairs.end())
      if(pairs.find(readname) != pairs.end())
	{
	      out << start << '\t'
		  << chrom << '\t' 
		  << strand << '\n';
	  ++nfound;
	  if(nfound == pairs.size()) return;
	}
    }
}

int main( int argc, char ** argv )
{
  int argn = 1;
  const char * umm_mapfile = argv[argn++]; // the file where all multis in U/M pairs map.  Part of output of process_unsorted_archive
  const char * umu_mapfile = argv[argn++]; // the file where all uniques in U/M pairs map.  Part of output of process_unsorted_archive
  const unsigned minmqual = atoi(argv[argn++]); //minimum mapping quality
  //read in where the TE are in the reference
  const char * outfile = argv[argn++];
  char * ref_data = NULL;
  if ( argc > 5 )
    {
      ref_data = argv[argn++]; // a tab-delimited file: chromosome start stop for each TE in reference
    }

  vector< teinfo > refdata;
  
  if( ref_data != NULL )
    {
      ifstream in(ref_data);
      unsigned chrom,start,stop;
      //string name;
      while ( !in.eof() )
	{
	  in >> chrom >> start >> stop >> ws;
	  refdata.push_back( teinfo(chrom,start,stop) );
	}
      in.close();
      cerr << "finished with reference data\n";
    }

  //set<unsigned> umm;
  filtering_ostream out;
  out.push(gzip_compressor());
  out.push(file_sink(outfile,ios_base::binary|ios_base::out));
  maptype umm;
  if( isbinary(umm_mapfile) )
    {
      filtering_istream gzin;
      gzin.push(gzip_decompressor());
      gzin.push(file_source(umm_mapfile));
      read_data(umm,
		refdata,
		gzin,minmqual);
    }
  else
    {
      ifstream in(umm_mapfile);
      read_data(umm,refdata,
		in,minmqual);
    }
  cerr << "finished um_m file\n";
  if ( isbinary(umu_mapfile) )
    {
      filtering_istream gzin;
      gzin.push(gzip_decompressor());
      gzin.push(file_source(umu_mapfile));
      read_umu(gzin,out,umm);
    }
  else
    {
      ifstream in(umu_mapfile);
      read_umu(in,out,umm);
    }
}
