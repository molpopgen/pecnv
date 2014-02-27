/*
  1. Assume bam file sorted by read name
  2. usage samtools -f 2 | bwa_mapdistance outfilename
  NOTE: outfile is in .gz format
*/

#include <Sequence/samrecord.hpp>
#include <map>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <bwa_util.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

using namespace std;
using namespace Sequence;
using namespace boost::iostreams;

int main( int argc, char ** argv )
{
  int argn = 1;
  const char * ofilename = argv[argn++];

  samrecord r1,r2;

  map<unsigned,unsigned> mdist;

  //unsigned lane1,rpair1,lane2,rpair2;
  unsigned nproc=0;
  while(!cin.eof())
    {
      cin >> r1 >> ws >> r2 >> ws;
      ++nproc;
      if(nproc % 1000000 == 0.)
	{
	  cerr << nproc << " processed\n";
	}
      /*
	getLanePair(&lane1,&rpair1,r1.qname());
	getLanePair(&lane2,&rpair2,r2.qname());
      */
      if ( r1.qname() != r2.qname() )
	{
	  cerr << "error: reads appear to not be properly sorted by read name:\n"
	       << r1 << '\n'
	       << r2 << '\n';
	  exit(1);
	}
      /*
      if( lane1 != lane2 || rpair1 != rpair2 )
	{
	  cerr << "error: lane/pair id mismatch:\n"
	       << r1 << '\n'
	       << r2 << '\n';
	  exit(1);
	}
      */
      if( r1.rname() == r2.rname() )
	{
	  //pair is linked
	  if(!hasXT(r1,"R") && !hasXT(r2,"R")) //if reads are unique and/or rescued by sampe
	    {
	      int mdist1 = r1.isize();
	      int pos1 = r1.pos(),pos2=r2.pos();
	      samflag f1 = r1.flag(),f2 = r2.flag();
	      if( 
		 ( pos1 < pos2 && ( (!f1.qstrand && f1.mstrand)
				    || ( f2.qstrand && !f2.mstrand) ) )
		 ||
		 ( ( pos2 < pos1 ) && ( (f1.qstrand && !f1.mstrand)
					|| ( !f2.qstrand && f2.mstrand ) ) )
		  )
		{
#ifndef NDEBUG
		  int mdist2=r2.isize();
		  assert( abs(mdist1) == abs(mdist2) );
#endif
		  map<unsigned,unsigned>::iterator itr =  mdist.find(abs(mdist1));
		  if( itr == mdist.end() )
		    {
		      mdist.insert(make_pair(abs(mdist1),1));
		    }
		  else
		    {
		      itr->second++;
		    }
		}
	      // else
	      // 	{
	      // 	  cerr << r1 << '\n'
	      // 	       << r2 << '\n';
	      // 	}
	    }
	}
    }
  unsigned sum = 0;
  for( map<unsigned,unsigned>::const_iterator i = mdist.begin(); 
       i != mdist.end() ; ++i )
    {
      sum += i->second;
    }
  filtering_ostream out;
  out.push(gzip_compressor());
  out.push(file_sink(ofilename),ios_base::binary|ios_base::out);
  out << "distance\tnumber\tcprob\n";
  unsigned cum=0;
  for( map<unsigned,unsigned>::const_iterator i = mdist.begin(); 
       i != mdist.end() ; ++i )
    {
      cum += i->second;
      out << i->first << '\t'
	  << i->second << '\t'
	  << scientific << double(cum)/double(sum) << '\n';
    }
  out.pop();
  out.pop();
}
