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
//using ReadCollection = vector< LeftRights >;
using ReadCollection = unordered_set<string>;
using PhrapInput = vector< pair<pair< vector<Fasta>, vector<Fasta> >,
				pair< vector<Fasta>, vector<Fasta> > > >;

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

  ReadCollection rv;//( cEs.size() );
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
	  if( cluster != cEs.cend() )
	    {
	      rv.insert(editRname(b.read_name()));
	    }
	}
    }
  
  return rv;
}

PhrapInput seqQual( const params & pars, const vector<clusteredEvent> & cEs,
		    const ReadCollection & r )
{
  bamreader reader(pars.bamfile.c_str());

  if(! reader )
    {
      cerr << "Error: could not open " << pars.bamfile
	   << " for reading. Line " << __LINE__
	   << " of " << __FILE__ << '\n';
      exit(1);
    }

  PhrapInput rv(cEs.size());
  while(!reader.eof() && !reader.error())
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break;
      string n = editRname(b.read_name());
      if ( r.find(n) != r.end() )
	{
	  samflag bf(b.flag()); 
	  int side = -1;//-1 = not known, 0 = left, 1 = right
	  string chrom = (reader.ref_cbegin() + b.refid())->first; //Get the string representing this read's chromosome
	  //Is the alignment start position of this read w/in ISIZE
	  //of a (filtered) cluster?
	  auto cluster  = find_if(cEs.cbegin(),cEs.cend(),
				  bind(Cfinder,placeholders::_1,bf.qstrand,b.pos(),pars.INSERTSIZE,chrom,&side));
	  while( cluster != cEs.cend() )
	    {
	      auto rvitr = rv.begin() + (cluster-cEs.cbegin());
	       	  string qstring;
		  for_each(b.qual_cbegin(),b.qual_cend(),
			   [&](const int8_t & __i) {
			     qstring += to_string(int(__i+33)) + ' ';
			   });
	      if(side==0)
		{
		  rvitr->first.first.push_back(Fasta(n,b.seq()));
		  rvitr->first.second.push_back(Fasta(n,qstring));
		}
	      else if (side == 1)
		{
		  rvitr->second.first.push_back(Fasta(n,b.seq()));
		  rvitr->second.second.push_back(Fasta(n,qstring));
		}
	      side = -1;
	      cluster = find_if(cluster+1,cEs.cend(),
				bind(Cfinder,placeholders::_1,bf.qstrand,b.pos(),pars.INSERTSIZE,chrom,&side));
	    }
	}
    }
  return rv;
}

string baseName(const string & basedir,
		const string & chrom,
		const int32_t & first,
		const int32_t & last,
		const short & side )
{
  string SIDE = (!side) ? ".left" : ".right";
  string n = basedir + "/" + chrom + "." + to_string(first) + "." + to_string(last) + SIDE + ".fasta";
  return n;
}

void write2file(const string & filename,
		const vector<Fasta> & vf )
{
  ofstream out(filename.c_str());
  if(!out) 
    {
      cerr << "Error: could not open "
	   << filename 
	   << " for writing at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }
  copy(vf.cbegin(),vf.cend(),ostream_iterator<const Fasta>(out,"\n"));
}

void output( const params & pars,
	     const vector<clusteredEvent> & cEs,
	     const PhrapInput & pI )
{
  auto __pi = pI.cbegin();
  for_each( cEs.cbegin(),
	    cEs.cend(),
	    [&](const clusteredEvent & __c) {
	      if( ! __pi->first.first.empty() ) //Then we have some "LEFT" reads stored
		{
		  assert( __c.pfirst != -1 && __c.plast != -1 );
		  assert( __pi->first.first.size() == __pi->first.second.size() );
		  string seqfilename = baseName(pars.phrapdir,__c.chrom,__c.pfirst,__c.plast,0);
		  string qualfilename = seqfilename + ".qual";
		  write2file(seqfilename,__pi->first.first);
		  write2file(qualfilename,__pi->first.second);
		}	      
	      if( ! __pi->second.first.empty() ) //Then we have some "LEFT" reads stored
		{
		  assert( __c.mfirst != -1 && __c.mlast != -1 );
		  assert( __pi->second.first.size() == __pi->second.second.size() );
		  string seqfilename = baseName(pars.phrapdir,__c.chrom,__c.mfirst,__c.mlast,1);
		  string qualfilename = seqfilename + ".qual";
		  write2file(seqfilename,__pi->second.first);
		  write2file(qualfilename,__pi->second.second);
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
  vector<clusteredEvent> cEs = parseClusters(clusters,pars);
 
  //Collect the names of all read pairs on the proper strand for each event
  auto LR = getRnames(pars,cEs);

  //Get the sequences and quality scores for all sequences for both sides
  auto sq = seqQual(pars,cEs,LR);

  //Output
  output(pars,cEs,sq);
}
