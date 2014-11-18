/*
  Estimates insert size distribution from "proper pairs"
  in a BAM file
*/

#include <Sequence/bamreader.hpp>
#include <Sequence/samflag.hpp>
#include <map>
#include <limits>
#include <unordered_map>
#include <iostream>
#include <zlib.h>

#include <common.hpp>
using namespace std;
using namespace Sequence;

int main( int argc, char ** argv )
{
  int argn = 1;
  if ( argc < 3 )
    {
      cerr << "Usage: "
	   << argv[0] << " bamfile output.gz MAXPAIRS\n"
	   << "Where MAXPAIRS = maximum number of read pairs to evaluate.\n"
	   << "MAXPAIRS is an optional argument.  The default is to use all alignments\n"
	   << "in the BAM file.\n";
      exit(0);
    }
  const char * bamfilename = argv[argn++];
  const char * ofilename = argv[argn++];

  unsigned MAXPAIRS = numeric_limits<unsigned>::max();
  if( argc == 4 )
    {
      MAXPAIRS = stoi(argv[argn]);
    }

  bamreader reader(bamfilename);
  if( ! reader )
    {
      cerr << "Error: " << bamfilename
	   << " could not be opened for reading.\n";
      exit(1);
    }

  map<unsigned,unsigned> mdist;
  unordered_map<string,bamrecord> reads;
  unsigned PAIRS_EVALUATED = 0;
  while(!reader.eof()&&!reader.error())
    {
      bamrecord b = reader.next_record();
      if(!b.empty())
	{
	  samflag sf = b.flag();
	  if( sf.is_proper_pair  )
	    {
	      bamaux a1 = b.aux("XT");
	      if(a1.size)
		{
		  if (a1.value[0] != 'R')
		    {
		      string n = editRname(b.read_name());
		      auto itr = reads.find(n);
		      if( itr == reads.end() )
			{
			  reads.insert(make_pair(move(n),move(b)));
			}
		      else
			{
			  if( b.refid() == itr->second.refid() )
			    {
			      samflag sf2 = itr->second.flag();
			      auto pos1 = b.pos(),pos2=itr->second.pos();
			      if(( pos1 < pos2 && ( (!sf.qstrand && sf.mstrand)
						    || ( sf2.qstrand && !sf2.mstrand) ) )
				 ||
				 ( ( pos2 < pos1 ) && ( (sf.qstrand && !sf.mstrand)
							|| ( !sf2.qstrand && sf2.mstrand ) ) )
				 )
				{
				  auto mitr = mdist.find(abs(b.tlen()));
				  if( mitr == mdist.end() )
				    {
				      mdist.insert(make_pair(abs(b.tlen()),1));
				    }
				  else
				    {
				      mitr->second++;
				      ++PAIRS_EVALUATED;
				    }
				}
			    }
			  reads.erase(itr);
			}
		    }
		}
	    }
	}
      if( MAXPAIRS != numeric_limits<unsigned>::max() && PAIRS_EVALUATED >= MAXPAIRS ) break;
    }

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
      if( gzprintf(out,"%u\t%u\t%e\n",i->first,i->second,double(cum)/double(sum)) <= 0 )
	{
	  cerr << "Error: gzprintf error at line "
	       << __LINE__ << " of " << __FILE__ << '\n';
	  exit(1);
	}
    }
  gzclose(out);
}
