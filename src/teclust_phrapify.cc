#include <teclust_phrapify.hpp>
#include <common.hpp>

#include <sstream>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <utility>
#include <functional>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <cassert>
#include <Sequence/bamreader.hpp>
#include <Sequence/samfunctions.hpp>
#include <Sequence/Fasta.hpp>

//For file creation, etc.
#include <sys/stat.h>
#include <cerrno>

using namespace std;
using namespace Sequence;

struct clusteredEvent
{
  string chrom;
  int32_t nplus,nminus,pfirst,plast,pdist,pin,mfirst,mlast,mdist,min;
  std::istream & read(istream &);
};

istream & clusteredEvent::read(istream & in)
{
  in >> chrom >> nplus >> nminus >> pfirst >> plast
     >> pdist >> pin >> mfirst >> mlast 
     >> mdist >> min;
  return in;
}

istream & operator>>(istream & in, clusteredEvent & ce)
{
  return ce.read(in);
}

vector<clusteredEvent> parseClusters(const string & clusters,const params & pars)
//Replaces functionality of filter_edit, but w/more flexibility
{
  vector<clusteredEvent> cEs;
  istringstream in(clusters);
  string temp;
  getline(in,temp);//get rid of header
  while(!in.eof()) //read records
    {
      clusteredEvent e;
      in >> e >> ws;
      if(e.nplus >= pars.MINREADS && e.nminus >= pars.MINREADS) //filter on read number
	{
	  if( e.pdist >= pars.CLOSEST && e.mdist >= pars.CLOSEST && 
	      (!pars.novelOnly || (pars.novelOnly && e.pin == 0 && e.min == 0)) )
	    {
	      cEs.push_back(e);
	    }
	}
    }
  return cEs;
}

using LeftRights = pair<unordered_set<string>,
			unordered_set<string> >;
using ReadCollection = vector< LeftRights >;
using PhrapInput = vector< pair<pair< vector<Fasta>, vector<Fasta> >,
				pair< vector<Fasta>, vector<Fasta> > > >;

ReadCollection
getRnames( const params & pars,
	   const vector<clusteredEvent> & cEs )
{

  bamreader reader(pars.bamfile.c_str());

  if(! reader )
    {
      cerr << "Error: could not open " << pars.bamfile
	   << " for reading. Line " << __LINE__
	   << " of " << __FILE__ << '\n';
      exit(1);
    }

  auto Cfinder = [](const clusteredEvent & __cE,
		    const bool & strand,
		    const int32_t & pos,
		    const unsigned & INSERTSIZE,
		    const string & chrom,
		    int * side)
    {
      if(__cE.chrom != chrom) return false;
      if( !strand ) //if read is on + strand
	{
	  //If read is 5' of this cluster's end
	  if( pos < __cE.plast &&
	      __cE.pfirst != -1 && __cE.plast != -1 &&
	      ( abs(pos-__cE.pfirst) <= INSERTSIZE ||
		abs(pos-__cE.plast) <= INSERTSIZE ) )
	    {
	      *side = 0; //LEFT
	      return true;
	    }
	}
      else if (strand) //read is on - strand
	{
	  //If read is 3' of this cluster's start
	  if( pos > __cE.mfirst &&
	      __cE.mfirst != -1 && __cE.mlast != -1 &&
	      ( abs(pos-__cE.mfirst) <= INSERTSIZE ||
		abs(pos-__cE.mlast) <= INSERTSIZE ) )
	    {
	      *side = 1; //RIGHT
	      return true;
	    }
	}
      return false;
    };

  ReadCollection rv( cEs.size() );
  while(!reader.eof() && !reader.error())
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break;
      samflag bf(b.flag());
      if(!bf.query_unmapped)
	{
	  int side = -1;//-1 = not known, 0 = left, 1 = right
	  string chrom = (reader.ref_cbegin() + b.refid())->first; //Get the string representing this read's chromosome
	  //Is the alignment start position of this read w/in ISIZE
	  //of a (filtered) cluster?
	  auto cluster  = find_if(cEs.cbegin(),cEs.cend(),
				  bind(Cfinder,placeholders::_1,bf.qstrand,b.pos(),pars.INSERTSIZE,chrom,&side));
	  while( cluster != cEs.cend() )
	    {
	      auto n = editRname(b.read_name());
	      auto rvitr = rv.begin() + (cluster-cEs.cbegin());
	      if( !side && rvitr->first.find(n) == rvitr->first.end())
		{
		  rvitr->first.insert(n);
		}
	      else if(side==1 && rvitr->second.find(n) == rvitr->second.end())
		{
		  rvitr->second.insert(n);
		}
	      side = -1;
	      cluster = find_if(cluster+1,cEs.cend(),
				 bind(Cfinder,placeholders::_1,bf.qstrand,b.pos(),pars.INSERTSIZE,chrom,&side));
	    }
	}
    }
  
  return rv;
}

PhrapInput seqQual( const params & pars, const ReadCollection & r )
{
  bamreader reader(pars.bamfile.c_str());

  if(! reader )
    {
      cerr << "Error: could not open " << pars.bamfile
	   << " for reading. Line " << __LINE__
	   << " of " << __FILE__ << '\n';
      exit(1);
    }

  auto FindCollection = [](const LeftRights & __lr, 
			   const string & n,
			   bool * isleft) 
    {
      if ( __lr.first.find(n) != __lr.first.cend() )
	{
	  return true;
	}
      else if ( __lr.second.find(n) != __lr.second.cend() )
	{
	  *isleft = false;
	  return true;
	}
      return false;
    };

  PhrapInput rv(r.size());
  while(!reader.eof() && !reader.error())
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break;
      samflag bf(b.flag());

      bool isleft = true;
      string n = editRname(b.read_name());
      auto itr = find_if( r.cbegin(),
			  r.cend(),
			  bind(FindCollection,placeholders::_1,n,&isleft) );
      //if( itr != r.cend() )
      while( itr != r.cend() )
	{
	  string qstring;
	  for_each(b.qual_cbeg(),b.qual_cend(),
		   [&](const int8_t & __i) {
		     qstring += (__i+33) + ' ';
		   });
	  auto rvItr = (rv.begin() + size_t(itr-r.cbegin()));
	  if (isleft)
	    {
	      rvItr->first.first.emplace_back(Fasta(n,b.seq()));
	      rvItr->first.second.emplace_back(Fasta(n,move(qstring)));
	    }
	  else
	    {
	      rvItr->second.first.emplace_back(Fasta(n,b.seq()));
	      rvItr->second.second.emplace_back(Fasta(n,move(qstring)));
	    }
	  isleft = true;
	  itr = find_if(itr+1,
			r.cend(),
			bind(FindCollection,placeholders::_1,n,&isleft) );
	}
    }
  return rv;
}

void output( const params & pars,
	     const vector<clusteredEvent> & cEs,
	     const PhrapInput & pI )
{
  auto __pi = pI.cbegin();
  for_each( cEs.cbegin(),
	    cEs.cbegin(),
	    [&](const clusteredEvent & __c) {
	      if( ! __pi->first.first.empty() ) //Then we have some "LEFT" reads stored
		{
		  assert( __c.pfirst != -1 && __c.plast != -1 );
		  assert( __pi->first.first.size() == __pi->first.second.size() );
		  string seqfilename = pars.phrapdir + "/" + __c.chrom + "." + to_string(__c.pfirst) + "." + to_string(__c.plast) + ".fasta";
		  string qualfilename = seqfilename + ".qual";

		  ofstream out(seqfilename.c_str());
		  if(!out) 
		    {
		      cerr << "Error: could not open "
			   << seqfilename 
			   << " for writing at line " << __LINE__ 
			   << " of " << __FILE__ << '\n';
		      exit(1);
		    }
		  //Write out the sequences
		  copy( __pi->first.first.cbegin(),
			__pi->first.first.cend(),
			ostream_iterator<const Fasta>(out,"\n") );
		  out.close();
		  out.open(qualfilename.c_str());
		  if(!out) 
		    {
		      cerr << "Error: could not open "
			   << seqfilename 
			   << " for writing at line " << __LINE__ 
			   << " of " << __FILE__ << '\n';
		      exit(1);
		    }
		  copy( __pi->first.second.cbegin(),
			__pi->first.second.cend(),
			ostream_iterator<const Fasta>(out,"\n") );
		  out.close();
		}	      
	      if( ! __pi->second.first.empty() ) //Then we have some "LEFT" reads stored
		{
		  assert( __c.mfirst != -1 && __c.mlast != -1 );
		  assert( __pi->second.first.size() == __pi->second.second.size() );
		  string seqfilename = pars.phrapdir + "/" + __c.chrom + "." + to_string(__c.mfirst) + "." + to_string(__c.mlast) + ".fasta";
		  string qualfilename = seqfilename + ".qual";
		  ofstream out(seqfilename.c_str());
		  if(!out) 
		    {
		      cerr << "Error: could not open "
			   << seqfilename 
			   << " for writing at line " << __LINE__ 
			   << " of " << __FILE__ << '\n';
		      exit(1);
		    }
		  //Write out the sequences
		  copy( __pi->second.first.cbegin(),
			__pi->second.first.cend(),
			ostream_iterator<const Fasta>(out,"\n") );
		  out.close();
		  out.open(qualfilename.c_str());
		  if(!out) 
		    {
		      cerr << "Error: could not open "
			   << seqfilename 
			   << " for writing at line " << __LINE__ 
			   << " of " << __FILE__ << '\n';
		      exit(1);
		    }
		  copy( __pi->second.second.cbegin(),
			__pi->second.second.cend(),
			ostream_iterator<const Fasta>(out,"\n") );
		  out.close();
		}
	      ++__pi;
	    });
}

void phrapify( const params & pars,
	       const string & clusters )
{
  if(pars.phrapdir.empty()) return;
  if(pars.bamfile.empty()) return;

  //Try to create the output directory
  int status = mkdir( pars.phrapdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
  if( status == -1 && errno != EEXIST )
    {
      cerr << "Error: could not create directory "
	   << pars.phrapdir 
	   << " at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }
  //This is the new "filter_edit"
  cerr << "Filtering...\n";
  vector<clusteredEvent> cEs = parseClusters(clusters,pars);
 
  //Collect the names of all read pairs on the proper strand for each event
  cerr << "Collecting...\n";
  auto LR = getRnames(pars,cEs);

  //Get the sequences and quality scores for all sequences for both sides
  cerr << "Retrieving sequences and quality scores...\n";
  auto sq = seqQual(pars,LR);

  //Output
  cerr << "Output...\n";
  output(pars,cEs,sq);
  cerr << "Done" << endl;
}
