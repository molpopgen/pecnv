/*
  1. Assume bam file sorted by read name
  2. usage samtools -f 2 | bwa_mapdistance outfilename
  NOTE: outfile is in .gz format
*/

#include <Sequence/bamreader.hpp>
#include <Sequence/samflag.hpp>
#include <map>
#include <unordered_map>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <zlib.h>
/*
#include <bwa_util.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
*/
using namespace std;
using namespace Sequence;
//using namespace boost::iostreams;

string editRname( const std::string & readname )
{
  auto pound = readname.find('#');
  if(pound == string::npos) return readname;
  return string(readname.begin(),readname.begin()+pound);
}

int main( int argc, char ** argv )
{
  int argn = 1;
  if ( argc != 3 )
    {
      cerr << "Usage: "
	   << argv[0] << " bamfile output.gz\n";
      exit(0);
    }
  const char * bamfilename = argv[argn++];
  const char * ofilename = argv[argn++];

  //unsigned lane1,rpair1,lane2,rpair2;
  unsigned nproc=0;
  bamreader reader(bamfilename);
  if( ! reader )
    {
      cerr << "Error: " << bamfilename
	   << " could not be opened for reading.\n";
      exit(1);
    }

  map<unsigned,unsigned> mdist;
  unordered_map<string,bamrecord> reads;
  while(!reader.eof()&&!reader.error())
    {
      bamrecord b = reader.next_record();
      if(!b.empty())
	{
	  samflag sf = b.flag();
	  auto pos1 = b.pos(),pos2=b.next_pos();
	  if( 
	     ( pos1 < pos2 && ( (!sf.qstrand && sf.mstrand)
				|| ( sf.qstrand && !sf.mstrand) ) )
	     ||
	     ( ( pos2 < pos1 ) && ( (sf.qstrand && !sf.mstrand)
				    || ( !sf.qstrand && sf.mstrand ) ) )
	      )
	    {
	      string n = editRname(b.read_name());
	      auto itr = reads.find(n);
	      if( itr == reads.end() )
		{
		  reads.insert(make_pair(move(n),move(b)));
		}
	      else
		{
		  //std::cerr << "found ";
		  bamaux a1 = b.aux("XT"),a2=itr->second.aux("XT");
		  if(a1.size && a2.size)
		    {
		      if(a1.value[0]!='R' && a2.value[0]!='R')
			{
			  //std::cerr << " made it\n";
			  auto mitr = mdist.find(abs(b.tlen()));
			  if(mitr == mdist.end())
			    {
			      mdist.insert(make_pair(abs(b.tlen()),1));
			    }
			  else 
			    mitr->second++;
			  // for_each( mdist.begin(),mdist.end(),
			  // 		[](const pair<unsigned,unsigned>&puu){
			  // 		  cerr << puu.first << '\t'<<puu.second << '\n';
			  // 		});
			}
		      // else 
		      //   {
		      //     cerr << "not right..." 
		      // 	   << a1.value[0] << ' ' << a2.value[0] << '\n';
		      //   }
		    }
		  reads.erase(itr);
		}
	    }
	}
    }

 //   while(!cin.eof())
 //     {
 //       cin >> r1 >> ws >> r2 >> ws;
 //       ++nproc;

 //       if ( r1.qname() != r2.qname() )
 // 	{
 // 	  cerr << "error: reads appear to not be properly sorted by read name:\n"
 // 	       << r1 << '\n'
 // 	       << r2 << '\n';
 // 	  exit(1);
 // 	}
 //       if( r1.rname() == r2.rname() )
 // 	{
 // 	  //pair is linked
// 	  if(!hasXT(r1,"R") && !hasXT(r2,"R")) //if reads are unique and/or rescued by sampe
// 	    {
// 	      int mdist1 = r1.isize();
// 	      int pos1 = r1.pos(),pos2=r2.pos();
// 	      samflag f1 = r1.flag(),f2 = r2.flag();
// 	      if( 
// 		 ( pos1 < pos2 && ( (!f1.qstrand && f1.mstrand)
// 				    || ( f2.qstrand && !f2.mstrand) ) )
// 		 ||
// 		 ( ( pos2 < pos1 ) && ( (f1.qstrand && !f1.mstrand)
// 					|| ( !f2.qstrand && f2.mstrand ) ) )
// 		  )
// 		{
// #ifndef NDEBUG
// 		  int mdist2=r2.isize();
// 		  assert( abs(mdist1) == abs(mdist2) );
// #endif
// 		  map<unsigned,unsigned>::iterator itr =  mdist.find(abs(mdist1));
// 		  if( itr == mdist.end() )
// 		    {
// 		      mdist.insert(make_pair(abs(mdist1),1));
// 		    }
// 		  else
// 		    {
// 		      itr->second++;
// 		    }
// 		}
// 	    }
// 	}
//     }
  unsigned sum = 0;
  for( map<unsigned,unsigned>::const_iterator i = mdist.begin(); 
       i != mdist.end() ; ++i )
    {
      sum += i->second;
    }
  gzFile out = gzopen(ofilename,"w");
  if(out==NULL)
    {
      cerr << "Error: could not open " << ofilename
	   << " for writing.\n";
      exit(1);
    }
  string header("distance\tnumber\tcprob\n");
  if(!gzwrite(out,header.c_str(),header.size()))
    {
      cerr << "Error: gzwrite error at line "
	   << __LINE__ << " of " << __FILE__ << '\n';
      exit(1);
    }
  unsigned cum=0;
  for( map<unsigned,unsigned>::const_iterator i = mdist.begin(); 
       i != mdist.end() ; ++i )
    {
      cum += i->second;
      string x;
      if( gzprintf(out,"%u\t%u\t%e\n",i->first,i->second,double(cum)/double(sum)) <= 0 )
	{
	  cerr << "Error: gzprintf error at line "
	       << __LINE__ << " of " << __FILE__ << '\n';
	  exit(1);
	}
      /*
      out << i->first << '\t'
  	  << i->second << '\t'
  	  << scientific << double(cum)/double(sum) << '\n';
      */
    }
  gzclose(out);
  // out.pop();
  // out.pop();
}
