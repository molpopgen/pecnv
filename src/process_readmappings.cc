/*
  Assumptions: bam file is sorted be read name and generted with BWA + samtools
  Also, this code has only been tested (and results independently validated with long-read sequencing)
  on alignments generated with bwa 0.5.9 using the command lines provided in Rogers et al.

  A merge of the features in bwa_bam_tomapfiles2.cc and bwa_mapdistance.cc
  1.  assumes alignment records are sorted by read name
  2.  assumes only paired reads are in the stream (e.g. samtools view -f 1 [bamfile] | ./this_program
*/
#include <Sequence/samrecord.hpp>
#include <Sequence/samfunctions.hpp>
#include <iostream>
#include <fstream>
#include <set>
#include <cstdlib>
#include <cmath>
#include <utility>
#include <map>
#include <vector>
#include <string>
#include <cctype>
#include <limits>
#include <vector>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <zlib.h>
#include <bwa_util.hpp>

using namespace std;
using namespace Sequence;

struct output_files
{
  enum MAPTYPE {DIV,PAR,UL,UMU,UMM};
  string structural_fn,structural_sam_fn,
    um_u_fn,um_m_fn,um_sam_fn;

  gzFile structural,um_u,um_m,
    structural_sam,um_sam;

  output_files(const char * structural_base, const char * um_base) :
    structural_fn(structural_base),
    structural_sam_fn(structural_base),
    um_u_fn(um_base),
    um_m_fn(um_base),
    um_sam_fn(um_base)
  {
    structural_fn += ".csv.gz";
    structural_sam_fn += ".sam.gz";
    um_u_fn += "_u.csv.gz";
    um_m_fn += "_m.csv.gz";
    um_sam_fn += ".sam.gz";

    structural = gzopen(structural_fn.c_str(),"w");
    if ( structural == NULL ) {
      cerr << "Error, could not open " << structural_fn
	   << " for writing\n";
      exit(1);
    }
    gzbuffer(structural,65536);

    structural_sam = gzopen(structural_sam_fn.c_str(),"w");
    if ( structural_sam == NULL ) {
      cerr << "Error, could not open " << structural_sam_fn
	   << " for writing\n";
      exit(1);
    }
    gzbuffer(structural_sam,65536);

    um_u = gzopen(um_u_fn.c_str(),"w");
    if ( um_u == NULL ) {
      cerr << "Error, could not open " << um_u_fn
	   << " for writing\n";
      exit(1);
    }
    gzbuffer(um_u,65536);

    um_m = gzopen(um_m_fn.c_str(),"w");
    if ( um_m == NULL ) {
      cerr << "Error, could not open " << um_m_fn
	   << " for writing\n";
      exit(1);
    }
    gzbuffer(um_m,65536);

    um_sam = gzopen(um_sam_fn.c_str(),"w");
    if ( um_sam == NULL ) {
      cerr << "Error, could not open " << um_sam_fn
	   << " for writing\n";
      exit(1);
    }
    gzbuffer(um_sam,65536);
  }

  ~output_files()
  {
    gzclose(structural);
    gzclose(structural_sam);
    gzclose(um_u);
    gzclose(um_m);
    gzclose(um_sam);
  }
  gzFile stream(const MAPTYPE & m)
  {
    if (m == UMU)
      {
	return um_u;
      }
    else if (m == UMM)
      {
	return um_m;
      }
    assert( m == DIV || m == PAR || m == UL );
    return structural;
  }
};

std::string mtype2string( const output_files::MAPTYPE & m )
{
  switch (m)
    {
    case output_files::DIV:
      return "DIV";
      break;
    case output_files::PAR:
      return "PAR";
      break;
    case output_files::UL:
      return "UL";
      break;
    case output_files::UMU:
      return "U";
      break;
    case output_files::UMM:
      return "M";
      break;
    };
  return "NA";
}

bool hasXA(const samrecord & r)
{
  for( samrecord::tag_iterator i = r.tag_begin();
       i!=r.tag_end();++i)
    {
      if(i->tag() == "XA")
	{
	  return true;
	}
    }
  return false;
}

vector< pair<char,
	     unsigned> > parse_cigar(const string & cigar);

unsigned mm(const unsigned & nm,
	    const vector< pair<char,
	    unsigned> > & cigar_data);
unsigned ngaps(const vector< pair<char,unsigned> > & cigar_data);
unsigned alen(const vector< pair<char,unsigned> > & cigar_data);
void outputM( gzFile out,
	      const samrecord & r );
void checkMap(const samrecord & r1,
	      const samrecord & r2,
	      const output_files::MAPTYPE & m,
	      output_files & of);

int main(int argc, char ** argv)
{
  int argn = 1;
  const char * structural_base = argv[argn++];
  const char * um_base = argv[argn++];
  const char * mdistfile = argv[argn++];
  struct output_files of(structural_base,um_base);
  
  samrecord r1,r2;
  map<unsigned,unsigned> mdist; //distribution of insert sizes
  while( !cin.eof())
    {
      cin >> r1 >> ws >> r2 >> ws;

      if ( r1.qname() != r2.qname() )
	{
	  cerr << "error: alignment records not properly sorted by read name\n"
	       << r1 << '\n'
	       << r2 << '\n';
	  exit(10);
	}
 
      samflag rf = r1.flag(),
	rf2 = r2.flag();
      if((!rf.query_unmapped && !rf.mate_unmapped) &&
	 (!rf2.query_unmapped && !rf2.mate_unmapped) ) //Fixed bug/issue #1 from git
	{
	  bool ppaired = false;
	  //check if reads in expected orientation.
	  if(r1.rname() == r2.rname() && rf.is_proper_pair)
	    {
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
		      ppaired = true;
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
		}
	    }
	  string qref(r1.rname()),mref(r1.mrnm());
	  if(mref != "=" && qref != mref) //UL
	    {
	      assert(!ppaired);
	      checkMap(r1,r2,output_files::UL,of);
	    }
	  else
	    {
	      if(r1.pos() != r1.mpos()) //don't map to same pos on reference
		{
		  if(rf.qstrand == rf.mstrand)
		    {
		      assert(!ppaired);
		      checkMap(r1,r2,output_files::PAR,of);
		    }
		  else if (rf.qstrand == 0 &&
			   r1.pos() > r1.mpos())
		    {
		      assert(!ppaired);
		      checkMap(r1,r2,output_files::DIV,of);
		    }
		  else if (rf.mstrand == 0 &&
			   r1.mpos() > r1.pos())
		    {
		      assert(!ppaired);
		      checkMap(r1,r2,output_files::DIV,of);
		    }
		  else if(!rf.is_proper_pair)//is a putative UM
		    {
		      if( !(hasXT(r1,"U") && hasXT(r2,"U")) )
			{
			  if( !(hasXT(r1,"R")&&  hasXT(r2,"R")) )
			    { 
			      checkMap(r1,r2,output_files::UMU,of);
			    }
			}
		    }
		  else if ( (hasXA(r1) && !hasXA(r2)) ||
			    (hasXA(r2) && !hasXA(r1)) )
		    {
		      //could be ppaired due to rescued read
		      //make sure both are not labelled XT:A:U or R
		      if( !(hasXT(r1,"U") && hasXT(r2,"U")) )
			{
			  if( !(hasXT(r1,"R")&&  hasXT(r2,"R")) )
			    { 
			      checkMap(r1,r2,output_files::UMU,of);
			    }
			}
		    }
		}
	    }
	}
    }
  //get sum of insert size dist
  long unsigned sum = 0;
  for( map<unsigned,unsigned>::const_iterator i = mdist.begin(); 
       i != mdist.end() ; ++i )
    {
      sum += unsigned(long(i->second));
    }

  gzFile mdistout = gzopen(mdistfile,"w");
  if( mdistout == NULL ) {
    cerr << "Error: could not open "
	 << mdistfile
	 << " for writing\n";
    exit(1);
  }
  if( gzprintf(mdistout,"%s\n","distance\tnumber\tcprob") <= 0 )
    {
      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }
  unsigned long cum=0;
  for( map<unsigned,unsigned>::const_iterator i = mdist.begin(); 
       i != mdist.end() ; ++i )
    {
      cum += unsigned(long(i->second));

      if( gzprintf(mdistout,"%u\t%u\t%e\n",i->first,i->second,double(cum)/double(sum)) <= 0 )
	{
	  cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
    }
  gzclose(mdistout);
}

void checkMap(const samrecord & r1,
	      const samrecord & r2,
	      const output_files::MAPTYPE & m,
	      output_files & of)
{
  //make sure no foul-ups get this far
  string name2 = r2.qname();
 
#ifndef NDEBUG
  assert( r1.qname() == r2.qname() );
#endif

  bool isU1 = hasXT(r1,"U"),
    isU2 = hasXT(r2,"U"),
    isR1 = hasXT(r1,"R"),
    isR2 = hasXT(r2,"R"),
    hasXA1 = hasXA(r1),
    hasXA2 = hasXA(r2);
  if(isR1 && isR2)
    {
      return;
    }
  bool written = false;
  if(isU1 && isU2 )
    {
      written = true;
      ostringstream obuffer;
      obuffer<< r1.qname() << '\t'
	     << r1.mapq() << '\t'
	     << r1.rname() << '\t'
	     << r1.pos()-1 << '\t'
	     << r1.pos() + alignment_length(r1) -1 << '\t'
	     << r1.flag().qstrand << '\t'
	     << mismatches(r1) << '\t'
	     << ngaps(r1) << '\t'
	     << mtype2string(m) << '\t'
	     << r2.mapq() << '\t'
	     << r2.rname() << '\t'
	     << r2.pos()-1 << '\t'
	     << r2.pos() + alignment_length(r2) -1 << '\t'
	     << r2.flag().qstrand << '\t'
	     << mismatches(r2) << '\t'
	     << ngaps(r2) << '\t'
	     << mtype2string(m);
      if ( gzprintf(of.stream(m),"%s\n",obuffer.str().c_str()) <= 0 )
	{
	  cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      obuffer.str(string());
      obuffer << r1 << '\t' << r2;
      if( gzprintf(of.structural_sam,"%s\n",obuffer.str().c_str()) <= 0)
	{
	  cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
    }
  else 
    {
      bool U1M2 = (( (isU1||isR1) && !hasXA1 ) && hasXA2);
      bool U2M1 = (( (isU2||isR2) && !hasXA2 ) && hasXA1);
      if(U1M2)
	{
	  written = true;
	  
	  ostringstream obuffer;
	  obuffer << r1.qname() << '\t'
		  << r1.mapq() << '\t'
		  << r1.rname() << '\t'
		  << r1.pos()-1 << '\t'
		  << r1.pos() + alignment_length(r1) - 2 << '\t'
		  << r1.flag().qstrand << '\t'
		  << mismatches(r1) << '\t'
		  << ngaps(r1);
	  if( gzprintf(of.stream(output_files::UMU),
		       "%s\n",
		       obuffer.str().c_str()) <= 0 )
	    {
	      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
		   << " of " << __FILE__ << '\n';
	      exit(1);
	    }
	  
	  outputM( of.stream(output_files::UMM),
		   r2 );

	  obuffer.str(string());
	  obuffer << r1 << '\n' << r2;
	  if( gzprintf(of.um_sam,"%s\n",obuffer.str().c_str()) <= 0 )
	    {
	      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
		   << " of " << __FILE__ << '\n';
	      exit(1);
	    }
	}
      else if(U2M1 )
	{
	  written = true;
	  outputM( of.stream(output_files::UMM),
		   r1 );
	  
	  ostringstream obuffer;
	  obuffer << r2.qname() << '\t'
		  << r2.mapq() << '\t'
		  << r2.rname() << '\t'
		  << r2.pos()-1 << '\t'
		  << r2.pos() + alignment_length(r2) - 2 << '\t'
		  << r2.flag().qstrand << '\t'
		  << mismatches(r2) << '\t'
		  << ngaps(r2);// << '\n';
	  if (gzprintf(of.stream(output_files::UMU),
		       "%s\n",
		       obuffer.str().c_str()) <= 0 )
	    {
	      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
		   << " of " << __FILE__ << '\n';
	      exit(1);
	    }
	  obuffer.str(string());
	  obuffer << r1 << '\n'<<r2;
	  if (gzprintf(of.um_sam,"%s\n",obuffer.str().c_str()) <= 0)
	    {
	      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
		   << " of " << __FILE__ << '\n';
	      exit(1);
	    }
	}
    }
}

unsigned alen(const vector< pair<char,
	      unsigned> > & cigar_data)
{
  unsigned sum=0;
  for(unsigned i=0;i<cigar_data.size();++i)
    {
      if ( cigar_data[i].first == 'M' ||
	   cigar_data[i].first == 'I' ||
	   cigar_data[i].first == 'D' ||
	   cigar_data[i].first == 'N' )
	{
	  sum += cigar_data[i].second;
	}
    }
  return sum;
}

unsigned mm(const unsigned & nm,
	    const vector< pair<char,
	    unsigned> > & cigar_data)
{
  return nm - ngaps(cigar_data);
}

unsigned ngaps(const vector< pair<char,
	       unsigned> > & cigar_data)
{
  unsigned sum=0;
  for(unsigned i=0;i<cigar_data.size();++i)
    {
      if ( cigar_data[i].first == 'I' ||
	   cigar_data[i].first == 'D' )
	{
	  sum += cigar_data[i].second;
	}
    }
  return sum;
}

vector<pair<char,
	    unsigned> > parse_cigar(const string & cigar)
{
  vector<pair<char,
	      unsigned> > cigar_data;
  
  string::const_iterator pibeg = find_if( cigar.begin(),cigar.end(),::isdigit);
  string::const_iterator piend = find_if( pibeg+1,cigar.end(),::isalpha );
  char * endptr;
  while( pibeg != cigar.end() )
    {
      cigar_data.push_back( make_pair( *piend, 
				       strtoul( string(pibeg,piend).c_str(),&endptr,10 ) ) );
      pibeg = find_if( piend+1, cigar.end(), ::isdigit );
      piend = find_if( pibeg+1,cigar.end(),::isalpha);
    }
  return cigar_data;
}

struct mapping_pos
{
  string chrom;
  unsigned start,stop,strand,mm,gap;
  mapping_pos(const string &__c,
	      const unsigned &__s,
	      const unsigned &__st,
	      const unsigned &__str,
	      const unsigned &__mm,
	      const unsigned &__gap) : chrom(__c),
				       start(__s),
				       stop(__st),
				       strand(__str),
				       mm(__mm),
				       gap(__gap)
  {
  }
};

bool operator==(const mapping_pos & left, 
		const mapping_pos & right)
{
  bool same =  (left.chrom == right.chrom &&
		left.start == right.start &&
		left.stop == right.stop &&
		left.strand == right.strand &&
		left.mm == right.mm &&
		left.gap == right.gap);
  return same;
}
bool operator<(const mapping_pos & left, 
	       const mapping_pos & right)
{
  return (left.chrom < right.chrom &&
	  left.start < right.start &&
	  left.strand < right.strand);
}

vector<mapping_pos> get_mapping_pos(const samrecord & r)
{
  vector<mapping_pos> rv;
  rv.push_back( mapping_pos( r.rname(), 
			     r.pos()-1,
			     r.pos()+alignment_length(r)-2,
			     r.flag().qstrand,
			     mismatches(r),
			     ngaps(r) ) );
  string XA;
  bool found=false;
  for( samrecord::tag_iterator i = r.tag_begin() ; 
       !found&&i != r.tag_end() ; ++i )
    {
      if(i->tag() == "XA")
	{
	  XA=i->value();
	  found=true;
	}
    }
  if ( found )
    {
      vector<string::size_type> colons;
      string::size_type colon = XA.find(";");
      do
	{
	  colons.push_back(colon);
	  colon = XA.find(";",colon+1);
	}
      while(colon != string::npos);
      string hit;
      for(unsigned i=0;i<colons.size();++i)
	{
	  if( i == 0 )
	    {
	      hit = string(XA.begin(),XA.begin()+colons[0]);
	    }
	  else
	    {
	      hit = string(XA.begin()+colons[i-1]+1,XA.begin()+colons[i]);
	    }
	  vector<string::size_type> commas;
	  string::size_type comma = hit.find(",");
	  do 
	    {
	      commas.push_back(comma);
	      comma = hit.find(",",comma+1);
	    }
	  while(comma != string::npos);

	  string hit_chrom = string(hit.begin(),hit.begin()+commas[0]);
	  int hit_start= atoi( string(hit.begin()+commas[0]+1,hit.begin()+commas[1]).c_str() );
	  string cigar(hit.begin()+commas[1]+1,hit.begin()+commas[2]);
	  vector<pair<char,unsigned> > cdata = parse_cigar(cigar);
	  unsigned hit_stop = abs(hit_start) + alen(cdata) - 2;
	  unsigned nm = atoi(string(hit.begin()+commas[2]+1,hit.end()).c_str());

	  mapping_pos hitmp( hit_chrom,abs(hit_start)-1,hit_stop, ((hit_start>0)?0:1),
			     mm(nm,cdata),ngaps(cdata));
	  if(find(rv.begin(),rv.end(),hitmp)==rv.end())
	    {
	      rv.push_back(hitmp);
	    }
	}
    }
  return rv;
}

void outputM( gzFile gzout,
	      const samrecord & r )
{
  string name = r.qname();
  vector<mapping_pos> mpos = get_mapping_pos(r);
  for( unsigned i=0;i<mpos.size();++i)
    {
      ostringstream out;
      out << name << '\t' << r.mapq() << '\t'
	  << mpos[i].chrom << '\t'
	  << mpos[i].start << '\t'
	  << mpos[i].stop << '\t'
	  << mpos[i].strand << '\t'
	  << mpos[i].mm << '\t'
	  << mpos[i].gap;// << '\n';
      if( gzprintf(gzout,"%s\n",out.str().c_str()) <= 0 )
	{
	  cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
    }
}
