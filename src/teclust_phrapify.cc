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
#include <cstring>
#include <map>
#include <thread>
#include <mutex>
#include <Sequence/bamreader.hpp>
#include <Sequence/samfunctions.hpp>
#include <Sequence/Fasta.hpp>

//For file creation, etc.
#include <sys/stat.h>
#include <cerrno>

using namespace std;
using namespace Sequence;

mutex output_mutex;

struct clusteredEvent
{
  string chrom;
  int32_t nplus,nminus,pfirst,plast,pdist,pin,mfirst,mlast,mdist,min;
  std::istream & read(istream &);     //Thru version 0.1.4
  std::istream & read_bed(istream &); //For BED format output, 0.1.5+
};

istream & clusteredEvent::read(istream & in)
{
  in >> chrom >> nplus >> nminus >> pfirst >> plast
     >> pdist >> pin >> mfirst >> mlast 
     >> mdist >> min;
  return in;
}

istream & clusteredEvent::read_bed(istream & in)
{
  string chr2,name,score,strand1,strand2,xtra;
  in >> chrom >> pfirst >> plast 
     >> chr2  >> mfirst >> mlast 
     >> name >> score >> strand1 >> strand2;
  getline(in,xtra);

  //Subtract 1 from plast and mlast to change to 0-offset, if they are > -1
  if( plast != -1 ) --pfirst;
  if( mlast != -1 ) --mlast;

  //if chrom is unknown (.), then chr2 must not be, so we can assign chr2 to chrom
  if(chrom == ".") chrom = chr2;
  assert( chrom != "." );

  //parse the extra stuff, which is 6 fields : nplus, nminus, pdist, pin, mdist, min
  char * tok = strtok( const_cast<char*>(xtra.c_str()), "\t" );
  short field = 0;
  while( tok != NULL )
    {
      if(field==0)
	{
	  nplus = stoi(tok);
	}
      else if(field==1)
	{
	  nminus = stoi(tok);
	}
      else if(field==2)
	{
	  pdist = stoi(tok);
	}
      else if(field==3)
	{
	  pin = stoi(tok);
	}
      else if(field==4)
	{
	  mdist = stoi(tok);
	}
      else if(field==5)
	{
	  min = stoi(tok);
	}
      ++field;
      tok = strtok(NULL,"\t");
    }
  if(field != 6)
    {
      cerr << "Error parsing xtra field data on line "
	   << __LINE__ << " of " << __FILE__ << '\n';
      exit(10);
    }
  return in;
}

istream & operator>>(istream & in, clusteredEvent & ce)
{
  //return ce.read(in);
  return ce.read_bed(in);
}

vector<clusteredEvent> parseClusters(const string & clusters,const teclust_params & pars)
//Replaces functionality of filter_edit, but w/more flexibility
{
  vector<clusteredEvent> cEs;
  istringstream in(clusters);
  string temp;
  //Assuming BEDPE input now, so no more header
  //getline(in,temp);//get rid of header
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
//using ReadCollection = unordered_set<string>;
using ReadCollection = unordered_set<int64_t>;
//using ReadCollection = unordered_set< pair<int32_t,int32_t> >;
using PhrapInput = vector< pair<pair< vector<Fasta>, vector<Fasta> >,
				pair< vector<Fasta>, vector<Fasta> > > >;

auto Cfinder = [](const clusteredEvent & __cE,
		  const bool & strand,
		  const int32_t & pos,
		  const int32_t & INSERTSIZE,
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
getRnames( const teclust_params & pars,
	   const vector<clusteredEvent> & cEs,
	   const int64_t * begin = nullptr,
	   const int64_t * end = nullptr)
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
  if( begin != nullptr ) reader.seek(*begin, SEEK_SET);
  while( (begin == nullptr && end == nullptr && !reader.eof() && !reader.error()) ||
	 (begin != nullptr && end != nullptr && reader.tell() <= *end) )
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
	      //rv.insert(editRname(b.read_name()));
	      int64_t c = b.refid();
	      c=(c<<32)|b.pos();
	      rv.insert(c);
	      //rv.insert(make_pair(b.refid(),b.pos()));
	    }
	}
    }
  
  return rv;
}

//for new threading approach
ReadCollection getRnames_v2( const bamrange * brange,
			     const teclust_params * pars,
			     const vector<clusteredEvent> * cEs )
{
  ReadCollection rv;

  bamreader reader(pars->bamfile.c_str());

  if(! reader )
    {
      cerr << "Error: could not open " << pars->bamfile
	   << " for reading. Line " << __LINE__
	   << " of " << __FILE__ << '\n';
      exit(1);
    }

  reader.seek(brange->beg,SEEK_SET);
  while( reader.tell() <= brange->end )
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break;
      auto bref = b.refid();
      if ( (b.refid()==brange->refid1 && b.pos() >= brange->start-1) ||
	   ( bref > brange->refid1 && bref < brange->refid2 ) ||
	   ( bref == brange->refid2 && b.pos() < brange->stop) )
	{
	  samflag bf(b.flag());
	  if(!bf.query_unmapped)
	    {
	      int side = -1;//-1 = not known, 0 = left, 1 = right
	      string chrom = (reader.ref_cbegin() + b.refid())->first; //Get the string representing this read's chromosome
	      //Is the alignment start position of this read w/in ISIZE
	      //of a (filtered) cluster?
	      auto cluster  = find_if(cEs->cbegin(),cEs->cend(),
				      bind(Cfinder,placeholders::_1,bf.qstrand,b.pos(),pars->INSERTSIZE,chrom,&side));
	      if( cluster != cEs->cend() )
		{
		  //rv.insert(editRname(b.read_name()));
		  int64_t c = b.refid();
		  c=(c<<32)|b.pos();
		  rv.insert(c);
		  //rv.insert(make_pair(b.refid(),b.pos()));
		}
	    }
	}
    }
  return rv;
  //Avoiding "false sharing"
  //rc = std::move(rv);
}

PhrapInput seqQual( const teclust_params & pars, const vector<clusteredEvent> & cEs,
		    const ReadCollection & r,
		    const int64_t * begin = nullptr,
		    const int64_t * end = nullptr)
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
  if(begin!=nullptr)reader.seek(*begin,SEEK_SET);
  //while(!reader.eof() && !reader.error())
  while( (begin == nullptr && end == nullptr && !reader.eof() && !reader.error()) ||
	 (begin != nullptr && end != nullptr && reader.tell() <= *end) )
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break;
      string n = editRname(b.read_name());
      //if ( r.find(n) != r.end() )
      int64_t c = b.refid();
      c=(c<<32);
      c|=b.pos();
      //if(r.find(make_pair(b.refid(),b.pos())) != r.end() )
      if(r.find(c) != r.end())
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

PhrapInput seqQual_v2( //PhrapInput & pi,
		 const bamrange * brange,
		 const teclust_params * pars, 
		 const vector<clusteredEvent> * cEs,
		 const ReadCollection * r )
{
   bamreader reader(pars->bamfile.c_str());

  if(! reader )
    {
      cerr << "Error: could not open " << pars->bamfile
	   << " for reading. Line " << __LINE__
	   << " of " << __FILE__ << '\n';
      exit(1);
    }

  PhrapInput rv(cEs->size());
  reader.seek(brange->beg,SEEK_SET);
  while( reader.tell() <= brange->end )
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break;
      auto bref = b.refid();
      if ( (b.refid()==brange->refid1 && b.pos() >= brange->start-1) ||
	   ( bref > brange->refid1 && bref < brange->refid2 ) ||
	   ( bref == brange->refid2 && b.pos() < brange->stop) )
	{
	  string n = editRname(b.read_name());
	  //if ( r.find(n) != r.end() )
	  int64_t c = b.refid();
	  c=(c<<32);
	  c|=b.pos();
	  //if(r.find(make_pair(b.refid(),b.pos())) != r.end() )
	  if(r->find(c) != r->end() )
	    {
	      samflag bf(b.flag()); 
	      int side = -1;//-1 = not known, 0 = left, 1 = right
	      string chrom = (reader.ref_cbegin() + b.refid())->first; //Get the string representing this read's chromosome
	      //Is the alignment start position of this read w/in ISIZE
	      //of a (filtered) cluster?
	      auto cluster  = find_if(cEs->cbegin(),cEs->cend(),
				      bind(Cfinder,placeholders::_1,bf.qstrand,b.pos(),pars->INSERTSIZE,chrom,&side));
	      while( cluster != cEs->cend() )
		{
		  auto rvitr = rv.begin() + (cluster-cEs->cbegin());
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
		  cluster = find_if(cluster+1,cEs->cend(),
				    bind(Cfinder,placeholders::_1,bf.qstrand,b.pos(),pars->INSERTSIZE,chrom,&side));
		}
	    }
	}
    }
  return rv;
  //pi = std::move(rv); 
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
  struct stat buf;
  bool exists = (stat(filename.c_str(),&buf) == 0) ? true:false;
  //open in append if file already exists
  ofstream out(filename.c_str(),(exists == true) ? ios_base::app : ios_base::out);
  //ofstream out(filename.c_str());//,(exists == true) ? ios_base::app : ios_base::out);
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

//I should just be able to lock with a mutex here, right??
void output( const teclust_params & pars,
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

void checkdir( const teclust_params & pars)
{
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
}

void phrapify( const teclust_params & pars,
	       const string & clusters )
{
  if(pars.phrapdir.empty()) return;
  if(pars.bamfile.empty()) return;

  checkdir(pars);

  //This is the new "filter_edit"
  vector<clusteredEvent> cEs = parseClusters(clusters,pars);
 
  //Collect the names of all read pairs on the proper strand for each event
  auto LR = getRnames(pars,cEs);

  //Get the sequences and quality scores for all sequences for both sides
  auto sq = seqQual(pars,cEs,LR);

  //Output
  output(pars,cEs,sq);
}

void phrapify_t_work(const teclust_params & pars, 
		     const vector<clusteredEvent> & cEs,
		     const int64_t & begin,
		     const int64_t & end)
{
  auto LR = getRnames(pars,cEs,&begin,&end);
  auto sq = seqQual(pars,cEs,LR,&begin,&end);
  output(pars,cEs,sq);
}

/*
  Threaded version.  1 thread per sequence in reference genome is the idea here
 */
void phrapify_t( const teclust_params & pars,
		 const vector<pair<string,pair<uint64_t,uint64_t> > > & offsets,
		 const string & clusters )
{
  if(offsets.empty())
    {
      cerr << "Warning from phrapify_t ( line " << __LINE__ 
	   << " of " __FILE__ << "): empty offsets vector encountered.  Defaulting to using a single thread";
      phrapify(pars,clusters);
    }
  if(pars.phrapdir.empty()) return;
  if(pars.bamfile.empty()) return;

  checkdir(pars);

  //This is the new "filter_edit"
  vector<clusteredEvent> cEs = parseClusters(clusters,pars);

  //rearrange cEs into a container that is per-reference arm.
  map<string,vector<clusteredEvent> > cEs2;
  for( auto i = cEs.begin() ; i != cEs.end() ; ++i )
    {
      cEs2[i->chrom].emplace_back(std::move(*i));
    }

  //vector<ReadCollection> vRC(cEs2.size());
  //vector<PhrapInput> vPI(cEs2.size());
  auto i = cEs2.begin();
  unsigned thread_id = 0;
  while( i != cEs2.end() )
    {
      vector<thread> vt(pars.NTHREADS);
      int32_t t = 0;
      for( ; t < pars.NTHREADS && i != cEs2.end(); ++t,++i )
	{
	  //find the chromo
	  auto __itr = find_if(offsets.begin(),offsets.end(),[&i](const pair<string,pair<uint64_t,uint64_t> > & __p) {
	      return i->first == __p.first;
	    });	  
	  vt[t]=thread(phrapify_t_work,cref(pars),cref(i->second),__itr->second.first,__itr->second.second);
	}
      for(int32_t t_i = 0 ; t_i < t ; ++t_i ) vt[t_i].join();
    }
}

/*
void getRnames_v2( ReadCollection & rc,
		   const bamrange & brange,
		   const teclust_params & pars,
		   const vector<clusteredEvent> & cEs )
void seqQual_v2( PhrapInput & pi,
		 const bamrange & brange,
		 const teclust_params & pars, const vector<clusteredEvent> & cEs,
		 const ReadCollection & r )
*/
void phrapify_t_work_v2(//ReadCollection & rc,
			//PhrapInput & pi,
			const bamrange * brange,
			const teclust_params * pars, 
			const vector<clusteredEvent> * cEs)
{
  ReadCollection rc = getRnames_v2(brange,pars,cEs);
  PhrapInput pi = seqQual_v2(brange,pars,cEs,&rc);
  /*
    prevent threads from writing to the same file,
    which could happen if the reference is divided up 
    into enough chunks
  */
  lock_guard<mutex> lockit(output_mutex);
  output(*pars,*cEs,pi);
}

void phrapify_t_v2( const teclust_params & pars,
		    const vector<bamrange> & branges,
		    const string & clusters )
{
  if (branges.empty())
    {
      cerr << "Warning from phrapify_t_v2 ( line " << __LINE__ 
	   << " of " __FILE__ << "): empty branges vector encountered.  Defaulting to using a single thread";
      phrapify(pars,clusters);
    }
 if(pars.phrapdir.empty()) return;
  if(pars.bamfile.empty()) return;

  checkdir(pars);

  //This is the new "filter_edit"
  vector<clusteredEvent> cEs = parseClusters(clusters,pars);
  vector<thread> vthreads(pars.NTHREADS);

  for( unsigned t = 0 ; t < vthreads.size() ; ++t )
    {
      vthreads[t] = thread(phrapify_t_work_v2,
			   &branges[t],&pars,&cEs);
    }
  for(unsigned t=0;t<vthreads.size();++t) vthreads[t].join();
}
