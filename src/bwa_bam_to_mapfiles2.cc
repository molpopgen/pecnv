/*
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
#include <cassert>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <bwa_util.hpp>

using namespace std;
using namespace boost::iostreams;
using namespace Sequence;

struct output_files
{
  enum MAPTYPE {DIV,PAR,UL,UMU,UMM};
  string structural_fn,structural_sam_fn,
    um_u_fn,um_m_fn,um_sam_fn;

  filtering_ostream structural,um_u,um_m,
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

    structural.push(gzip_compressor());
    structural.push(file_sink(structural_fn.c_str()),ios_base::out|ios_base::binary);

    structural_sam.push(gzip_compressor());
    structural_sam.push(file_sink(structural_sam_fn.c_str()),ios_base::out|ios_base::binary);

    um_u.push(gzip_compressor());
    um_u.push(file_sink(um_u_fn.c_str()),ios_base::out|ios_base::binary);

    um_m.push(gzip_compressor());
    um_m.push(file_sink(um_m_fn.c_str()),ios_base::out|ios_base::binary);

    um_sam.push(gzip_compressor());
    um_sam.push(file_sink(um_sam_fn.c_str()),ios_base::out|ios_base::binary);
  }

  ~output_files()
  {
    structural.pop();
    structural.pop();
    structural_sam.pop();
    structural_sam.pop();
    um_u.pop();
    um_u.pop();
    um_m.pop();
    um_m.pop();
    um_sam.pop();
    um_sam.pop();
  }
  filtering_ostream & stream(const MAPTYPE & m)//,
  //const string & readname)
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
    /*
    if(readname[readname.length()-1]=='0')
      {
	return structural_left;
      }
    return structural_right;
    */
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

/*
string fix_readname( const std::string & n )
{
  string rv(n);
  string::size_type pos = rv.find(':');
  while(pos != string::npos)
    {
      rv[pos]='\t';
      pos = rv.find(':',pos+1);
    }
  return rv;
}
*/
vector< pair<char,
	     unsigned> > parse_cigar(const string & cigar);

unsigned mm(const unsigned & nm,
	    const vector< pair<char,
			       unsigned> > & cigar_data);
unsigned ngaps(const vector< pair<char,
				  unsigned> > & cigar_data);
unsigned alen(const vector< pair<char,
				 unsigned> > & cigar_data);
filtering_ostream & outputM( filtering_ostream & out,
			     const samrecord & r );

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
  /*
    #ifndef NDEBUG
    unsigned lane1,rpair1,lane2,rpair2;
    getLanePair(&lane1,&rpair1,r1.qname());
    getLanePair(&lane2,&rpair2,r2.qname());
    assert(lane1==lane2&&rpair1==rpair2);
    #endif
  */
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
      of.stream(m)//,r1.qname()) 
	<< r1.qname() << '\t'
	//<< fix_readname(r1.qname()) << '\t' 
	<< r1.mapq() << '\t'
	<< r1.rname() << '\t'
	<< r1.pos()-1 << '\t'
	<< r1.pos() + alignment_length(r1) -1 << '\t'
	<< r1.flag().qstrand << '\t'
	<< mismatches(r1) << '\t'
	<< ngaps(r1) << '\t'
	<< mtype2string(m) << '\n';

      of.stream(m)//,r2.qname()) 
	<< r2.qname() << '\t'
	//<< fix_readname(r2.qname()) << '\t' 
	<< r2.mapq() << '\t'
	<< r2.rname() << '\t'
	<< r2.pos()-1 << '\t'
	<< r2.pos() + alignment_length(r2) -1 << '\t'
	<< r2.flag().qstrand << '\t'
	<< mismatches(r2) << '\t'
	<< ngaps(r2) << '\t'
	<< mtype2string(m) << '\n';

      of.structural_sam << r1 << '\n'
			<< r2 << '\n';
    }
  else //if((isU1||isR1) && isMorR2)//1 is R, 2 is U
    {
      bool U1M2 = (( (isU1||isR1) && !hasXA1 ) && hasXA2);
      bool U2M1 = (( (isU2||isR2) && !hasXA2 ) && hasXA1);
      if(U1M2)
	{
	  written = true;
	  
	  of.stream(output_files::UMU)//,r1.qname()) 
	    << r1.qname() << '\t'
	    //<< fix_readname(r1.qname()) << '\t' 
	    << r1.mapq() << '\t'
	    << r1.rname() << '\t'
	    << r1.pos()-1 << '\t'
	    << r1.pos() + alignment_length(r1) - 2 << '\t'
	    << r1.flag().qstrand << '\t'
	    << mismatches(r1) << '\t'
	    << ngaps(r1) << '\n';
	  
	  outputM( of.stream(output_files::UMM),//r2.qname()),
		   r2 );
	  of.um_sam << r1 << '\n'
		    << r2 << '\n';
	}
      //else if((isU2||isR2) && isMorR1)
      else if(U2M1 )
	{
	  //cerr << r2 << '\n' << r2.pos() <<'\n';
	  //cerr << r1.pos() << '\n';
	  written = true;
	  outputM( of.stream(output_files::UMM),//r1.qname()),
		   r1 );
	  
	  of.stream(output_files::UMU)//,r2.qname()) 
	    << r2.qname() << '\t'
	    //<< fix_readname(r2.qname()) << '\t' 
	    << r2.mapq() << '\t'
	    << r2.rname() << '\t'
	    << r2.pos()-1 << '\t'
	    << r2.pos() + alignment_length(r2) - 2 << '\t'
	    << r2.flag().qstrand << '\t'
	    << mismatches(r2) << '\t'
	    << ngaps(r2) << '\n';
	  
	  of.um_sam << r1 << '\n'
		    << r2 << '\n';
	}
    }
}


int main(int argc, char ** argv)
{
  int argn = 1;
  const char * structural_base = argv[argn++];
  const char * um_base = argv[argn++];
  struct output_files of(structural_base,um_base);
  
  samrecord r1,r2;
  //unsigned lane1,rpair1,lane2,rpair2;
  while( !cin.eof())
    {
      cin >> r1 >> ws >> r2 >> ws;
      //modification to move away from read name changing
      if ( r1.qname() != r2.qname() )
	{
	  cerr << "error: alignment records not properly sorted by read name\n"
	       << r1 << '\n'
	       << r2 << '\n';
	  exit(10);
	}
      /*
      getLanePair(&lane1,&rpair1,r1.qname());
      getLanePair(&lane2,&rpair2,r2.qname());
      if ( lane1 != lane2 ||
	   rpair1 != rpair2 )
	{
	  cerr << "error: alignment records not properly sorted by read name\n"
	       << r1 << '\n'
	       << r2 << '\n';
	  exit(10);
	}
      samflag rf = r1.flag(),
	rf2 = r2.flag();
      //Check BOTH flags.  Addresses isue #1 on github
      if( (!rf.query_unmapped && !rf.mate_unmapped) &&
	  (!rf2.query_unmapped && !rf2.mate_unmapped) )
	{
	  string qref(r1.rname()),mref(r1.mrnm());
	  if(mref != "=" && qref != mref) //UL
	    {
	      checkMap(r1,r2,output_files::UL,of);
	    }
	  else
	    {
	      if(r1.pos() != r1.mpos()) //don't map to same pos on reference
		{
		  if(rf.qstrand == rf.mstrand)
		    {
		      checkMap(r1,r2,output_files::PAR,of);
		    }
		  else if (rf.qstrand == 0 &&
			   r1.pos() > r1.mpos())
		    {
		      checkMap(r1,r2,output_files::DIV,of);
		    }
		  else if (rf.mstrand == 0 &&
			   r1.mpos() > r1.pos())
		    {
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
}

string::size_type parsingXA( const string & XA,
			     string::size_type & p,
			     string & chrom,
			     int & pos,
			     string & cigar,
			     unsigned & nm)
{
  if(p>XA.size()) return string::npos;

  string::size_type colon = XA.find(";",p);
  if(colon==string::npos) { return colon; }

  string::size_type comma = XA.find(",",p);
  assert(comma != string::npos);
  //chrom = atoi(string(XA.begin()+p,XA.begin()+comma).c_str());
  chrom = string(XA.begin()+p,XA.begin()+comma);
  p=comma+1;
  comma = XA.find(",",p);
  char sign = *(XA.begin()+p);
  assert(comma != string::npos);
  pos = (sign == '+') ? 
    atoi(string(XA.begin()+p+1,XA.begin()+comma).c_str()) :
    -atoi(string(XA.begin()+p+1,XA.begin()+comma).c_str());
  p=comma+1;
  comma = XA.find(",",p);
  assert(comma != string::npos);
  cigar = string(XA.begin()+p,XA.begin()+comma);
  p=comma+1;
  nm = atoi(string(XA.begin()+p,XA.begin()+colon).c_str());
  return ( (colon+1) < XA.length() ? (colon+1) : string::npos );
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
/*
filtering_ostream & outputM( filtering_ostream & out,
			     const samrecord & r )
{
  string name = fix_readname(r.qname());
  out << name << '\t'
      << r.mapq() << '\t'
      << r.rname() << '\t'
      << r.pos()-1 << '\t'
      << r.pos() + alignment_length(r) - 2 << '\t'
      << r.flag().qstrand << '\t'
      << mismatches(r) << '\t'
      << ngaps(r) << '\n';

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

  unsigned chrom,nm;
  int pos;
  string cigar;
  string::size_type p=0;
  while( (p=parsingXA(XA,p,chrom,pos,cigar,nm))!=string::npos )
    {
      vector<pair<char,unsigned> > cdata = parse_cigar(cigar);
      out << name << '\t'
	  << r.mapq() << '\t'
	  << chrom << '\t'
	//<< r.rname() << '\t'
	  << abs(pos)-1 << '\t'
	  << abs(pos) + alen(cdata) - 2 << '\t'
	  << ( (pos > 0) ? 0 : 1 ) << '\t'
	  << mm(nm,cdata) << '\t'
	  << ngaps(cdata) << '\n';
    }
  return out;
}
*/

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
  rv.push_back( mapping_pos( r.rname(), //atoi(r.rname().c_str()),
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
      //cerr << r << '\n';
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
	  //unsigned hit_chrom = atoi( string(hit.begin(),hit.begin()+commas[0]).c_str() );
	  string hit_chrom = string(hit.begin(),hit.begin()+commas[0]);
	  int hit_start= atoi( string(hit.begin()+commas[0]+1,hit.begin()+commas[1]).c_str() );
	  string cigar(hit.begin()+commas[1]+1,hit.begin()+commas[2]);
	  vector<pair<char,unsigned> > cdata = parse_cigar(cigar);
	  unsigned hit_stop = abs(hit_start) + alen(cdata) - 2;
	  unsigned nm = atoi(string(hit.begin()+commas[2]+1,hit.end()).c_str());
	  /*cerr << hit_chrom << ' '
	       << abs(hit_start) << ' ' 
	       << hit_stop << ' ' 
	       << nm << '\n';
	  */
	  mapping_pos hitmp( hit_chrom,abs(hit_start)-1,hit_stop, ((hit_start>0)?0:1),
			     mm(nm,cdata),ngaps(cdata));
	  if(find(rv.begin(),rv.end(),hitmp)==rv.end())
	    {
	      rv.push_back(hitmp);
	    }
	}
    }
  //cerr << "size = "<<rv.size()<<'\n';
  return rv;
}

filtering_ostream & outputM( filtering_ostream & out,
			     const samrecord & r )
{
  string name = r.qname();//fix_readname(r.qname());
  vector<mapping_pos> mpos = get_mapping_pos(r);
  //cerr << "name = " << name <<' ' << r.mapq() << '\n';
  for( unsigned i=0;i<mpos.size();++i)
    {
      out << name << '\t' << r.mapq() << '\t'
	  << mpos[i].chrom << '\t'
	  << mpos[i].start << '\t'
	  << mpos[i].stop << '\t'
	  << mpos[i].strand << '\t'
	  << mpos[i].mm << '\t'
	  << mpos[i].gap << '\n';
    }
  return out;
  /*
  out << name << '\t'
      << r.mapq() << '\t'
      << r.rname() << '\t'
      << r.pos()-1 << '\t'
      << r.pos() + alignment_length(r) - 2 << '\t'
      << r.flag().qstrand << '\t'
      << mismatches(r) << '\t'
      << ngaps(r) << '\n';

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

  unsigned chrom,nm;
  int pos;
  string cigar;
  string::size_type p=0;
  while( (p=parsingXA(XA,p,chrom,pos,cigar,nm))!=string::npos )
    {
      vector<pair<char,unsigned> > cdata = parse_cigar(cigar);
      out << name << '\t'
	  << r.mapq() << '\t'
	  << chrom << '\t'
	//<< r.rname() << '\t'
	  << abs(pos)-1 << '\t'
	  << abs(pos) + alen(cdata) - 2 << '\t'
	  << ( (pos > 0) ? 0 : 1 ) << '\t'
	  << mm(nm,cdata) << '\t'
	  << ngaps(cdata) << '\n';
    }
  return out;
  */
}
