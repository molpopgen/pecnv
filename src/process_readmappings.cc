/*
  Assumptions: bam file is sorted be read name and generted with BWA + samtools
  Also, this code has only been tested (and results independently validated with long-read sequencing)
  on alignments generated with bwa 0.5.9 using the command lines provided in Rogers et al.

  A merge of the features in bwa_bam_tomapfiles2.cc and bwa_mapdistance.cc
  1.  assumes alignment records are sorted by read name
  2.  assumes only paired reads are in the stream (e.g. samtools view -f 1 [bamfile] | ./this_program
*/
#include <Sequence/bamreader.hpp>
#include <Sequence/bamrecord.hpp>
#include <Sequence/samfunctions.hpp>
#include <iostream>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <zlib.h>

using namespace std;
using namespace Sequence;

using APAIR = pair<bamrecord,bamrecord>;
using readbucket = unordered_map<string, bamrecord>; //name, alignment
using PPairData = unordered_map<string,pair<bool,std::int32_t> >; //name, isU, tlen

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

    structural_sam = gzopen(structural_sam_fn.c_str(),"w");
    if ( structural_sam == NULL ) {
      cerr << "Error, could not open " << structural_sam_fn
	   << " for writing\n";
      exit(1);
    }

    um_u = gzopen(um_u_fn.c_str(),"w");
    if ( um_u == NULL ) {
      cerr << "Error, could not open " << um_u_fn
	   << " for writing\n";
      exit(1);
    }

    um_m = gzopen(um_m_fn.c_str(),"w");
    if ( um_m == NULL ) {
      cerr << "Error, could not open " << um_m_fn
	   << " for writing\n";
      exit(1);
    }

    um_sam = gzopen(um_sam_fn.c_str(),"w");
    if ( um_sam == NULL ) {
      cerr << "Error, could not open " << um_sam_fn
	   << " for writing\n";
      exit(1);
    }
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

vector< pair<char,unsigned> > parse_cigar(const string & cigar);

//No. mismatches
unsigned mm(const unsigned & nm,
	    const vector< pair<char,unsigned> > & cigar_data);
//No. gaps
unsigned ngaps(const vector< pair<char,unsigned> > & cigar_data);
//Alignment length
unsigned alen(const vector< pair<char,unsigned> > & cigar_data);
//Write U reads in U/P pair to files
void outputU( gzFile gzout,
	      gzFile gzoutSAM,
	      const bamrecord & r,
	      const bamreader & reader );
//Write M reads in U/P pair to files
void outputM( gzFile out,
	      gzFile gzoutSAM,
	      const bamrecord & r,
	      const bamreader & reader);
//Turn a read name from XXX#something to XXX
string editRname( const std::string & readname );



/*
  Does this pair of alignments represent a unique/multi pair?
*/
void evalUM(const bamrecord & b1,
	    const bamrecord & b2,
	    const bamreader & reader,
	    gzFile uout, gzFile mout,
	    gzFile SAMout);

void updateBucket( readbucket & rb, bamrecord & b, 
		   gzFile csvfile, gzFile samfile,
		   const char * maptype,
		   const bamreader & reader );

string toSAM(const bamrecord & b,
	     const bamreader & reader);

int main(int argc, char ** argv)
{
  int argn = 1;
  const char * bamfile = argv[argn++];
  const char * structural_base = argv[argn++];
  const char * um_base = argv[argn++];
  const char * mdistfile = argv[argn++];
  struct output_files of(structural_base,um_base);
  
  bamreader reader(bamfile);

  if ( ! reader ) {
    cerr << "Error: " << bamfile 
	 << " could not be opened for reading\n";
    exit(1);
  }

  readbucket DIV,PAR,UL,UM;
  PPairData ISIZES;
  map<unsigned,unsigned> mdist; //distribution of insert sizes
  auto pos = reader.tell(); //After the headers, @ start of 1st alignment
  while( !reader.eof() && !reader.error() )
    {
      bamrecord b = reader.next_record();
      samflag sf(b.flag());
      if(!sf.query_unmapped  && !sf.mate_unmapped) //Both reads are mapped
      	{
	  bamaux bXT = b.aux("XT");  //look for XT tag
	  if(bXT.size) //if it is present
	    {
	      const char XTval = bXT.value[0];
	      bool unusual = false; //A putative DIV/PAR/UL?
	      //Look for unusual read mappings here
	      if(XTval == 'U') //if read is uniquely-mapping
		{
		  if( b.refid() != b.next_refid() ) //both map to different scaffolds
		    {
		      unusual = true;
		      updateBucket(UL,b,of.structural,of.structural_sam,"UL\0",reader);
		    }
		  else if ( b.pos() != b.next_pos()) //Don't map to same position
		    {
		      if( sf.qstrand == sf.mstrand )
			{
			  unusual = true;
			  updateBucket(PAR,b,of.structural,of.structural_sam,"PAR\0",reader);
			}
		      else if( (sf.qstrand == 0 && b.pos() > b.next_pos()) ||
			       (sf.mstrand == 0 && b.next_pos() > b.pos() ) )
			{
			  unusual = true;
			  updateBucket(DIV,b,of.structural,of.structural_sam,"DIV\0",reader);
			}
		    }


		  if(!unusual)
		    {
		       string n = editRname(b.read_name());
		       auto i = UM.find(n);
		       if( i != UM.end() )
			 {
			   //Let's process the M/U pair and then delete it
			   //b is the unique-read, and the read
			   //at position i->second is the M/R read
			   //bamrecord multi = reader.record_at_pos(i->second);
			   //assert(!multi.empty());
			   evalUM(b,i->second,reader,of.um_u,of.um_m,of.um_sam);
			   UM.erase(i);
			 }
		       else 
			 {
			 }
		    }
		}
	      //putative U/M pair member, reads don't hit same position on same chromo
	      else if ((XTval == 'R' || XTval == 'M') && 
		       ( (b.refid() != b.next_refid()) ||
			 (b.refid() == b.next_refid() && b.pos() != b.next_pos()) ) )
		{
		  string n = editRname(b.read_name());
		  auto i = UM.find(editRname(n));
		  if(i == UM.end())
		    {				
		      UM.insert(make_pair(move(n),move(b)));
		    }
		  else //This is an M/M or M/R pair, so we can evaluate and then delete
		    {
		      evalUM(b,i->second,reader,of.um_u,of.um_m,of.um_sam);
		      UM.erase(i);
		    }
		}
	    }
	}
    }

  //These are done.
  DIV.clear();
  PAR.clear();
  UL.clear();

  if(!UM.empty())
    {
      cerr << "First pass complete, rewinding to finish scanning for U/M pairs\n";
      reader.seek( pos, SEEK_SET );
      
      while(!reader.eof() && !reader.error()) //This may not be working
	{
	  bamrecord b(reader.next_record());
	  if(b.empty()) break;
	  samflag r(b.flag());
	  if(!r.query_unmapped)
	    {
	      bamaux ba = b.aux("XT");
	      if(ba.value[0]=='U' || ba.value[0]=='R') //Read is flagged as uniquely-mapping or rescued
		{
		  string n = editRname(b.read_name());
		  auto i = UM.find(n);
		  if(i != UM.end()) //then the Unique reads redundant mate exists
		    {
		      evalUM(b,i->second,reader,of.um_u,of.um_m,of.um_sam);
		    }
		}
	    }
	}
    }
   
  //get sum of insert size dist
  // long unsigned sum = 0;
  // for( map<unsigned,unsigned>::const_iterator i = mdist.begin(); 
  //      i != mdist.end() ; ++i )
  //   {
  //     sum += unsigned(long(i->second));
  //   }

  // gzFile mdistout = gzopen(mdistfile,"w");
  // if( mdistout == NULL ) {
  //   cerr << "Error: could not open "
  // 	 << mdistfile
  // 	 << " for writing\n";
  //   exit(1);
  // }
  // if( gzprintf(mdistout,"%s\n","distance\tnumber\tcprob") <= 0 )
  //   {
  //     cerr << "Error: gzprintf error encountered at line " << __LINE__ 
  // 	   << " of " << __FILE__ << '\n';
  //     exit(1);
  //   }
  // unsigned long cum=0;
  // for( map<unsigned,unsigned>::const_iterator i = mdist.begin(); 
  //      i != mdist.end() ; ++i )
  //   {
  //     cum += unsigned(long(i->second));

  //     if( gzprintf(mdistout,"%u\t%u\t%e\n",i->first,i->second,double(cum)/double(sum)) <= 0 )
  // 	{
  // 	  cerr << "Error: gzprintf error encountered at line " << __LINE__ 
  // 	       << " of " << __FILE__ << '\n';
  // 	  exit(1);
  // 	}
  //   }
  // gzclose(mdistout);
}


string editRname( const std::string & readname )
{
  auto pound = readname.find('#');
  if(pound == string::npos) return readname;
  return string(readname.begin(),readname.begin()+pound);
}

void evalUM(const bamrecord & b1,
	    const bamrecord & b2,
	    const bamreader & reader,
	    gzFile uout, gzFile mout,
	    gzFile SAMout)
{
  if( editRname(b1.read_name()) != editRname(b2.read_name()) )
    {
      cerr << "Error: read names don't match at line "
	   << __LINE__
	   << " of " << __FILE__ << 'n';
      exit(1);
    }

  bamaux XTb1 = b1.aux("XT"),XTb2 = b2.aux("XT");
  if( XTb1.size && XTb2.size )
    {
      const char XTv1 = XTb1.value[0],
	XTv2 = XTb2.value[0];
      if ( (XTv1 == 'M' && XTv2 == 'M') ||
	   (XTv1 == 'R' && XTv2 == 'R') ) 
	return;
      bool U1M2 = ( ((XTv1=='U'||XTv1=='R') && b1.hasTag("XA")==nullptr) && b2.hasTag("XA") != nullptr );
      bool U2M1 = ( ((XTv2=='U'||XTv2=='R') && b2.hasTag("XA")==nullptr) && b1.hasTag("XA") != nullptr );
      if(U1M2)
	{
	  outputU(uout,SAMout,b1,reader);
	  outputM(mout,SAMout,b2,reader);
	  assert( !(XTv1=='M' && XTv2 == 'M') );
	}
      else if (U2M1)
	{
	  outputU(uout,SAMout,b2,reader);
	  outputM(mout,SAMout,b1,reader);
	  assert( !(XTv1=='M' && XTv2 == 'M') );
	}
    }
}

unsigned alen(const vector< pair<char,unsigned> > & cigar_data)
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

unsigned mm(const unsigned & nm,  const vector< pair<char, unsigned> > & cigar_data)
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

vector<pair<char,unsigned> > parse_cigar(const string & cigar)
{
  vector<pair<char,unsigned> > cigar_data;
  
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


//FXN NEED AUDITING
//The XA positions need to be turned into 0 offset
vector<mapping_pos> get_mapping_pos(const bamrecord & r,
				    const bamreader & reader)
{
  vector<mapping_pos> rv;
  auto REF = reader.ref_cbegin() + r.refid();
  rv.push_back( mapping_pos( REF->first,
			     r.pos(),
			     r.pos()+alignment_length(r)-1,
			     r.flag().qstrand,
			     mismatches(r),
			     ngaps(r) ) );
  bamaux auxXA = r.aux("XA");
  if(auxXA.size)
    {
      string XA(auxXA.value);
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

void outputU( gzFile gzout,
	      gzFile gzoutSAM,
	      const bamrecord & r,
	      const bamreader & reader )
{
  assert( ! r.flag().query_unmapped );
  assert( ! r.flag().mate_unmapped );
  string name = editRname(r.read_name());
  ostringstream obuffer;
  auto REF = reader.ref_cbegin()+r.refid();
  obuffer << name << '\t'
	  << r.mapq() << '\t'
	  << REF->first << '\t'
	  << r.pos() << '\t'
	  << r.pos() + alignment_length(r) - 1 << '\t'
	  << r.flag().qstrand << '\t'
	  << mismatches(r) << '\t'
	  << ngaps(r) << '\n'; 

  if(!gzwrite(gzout,obuffer.str().c_str(),obuffer.str().size()))
    {
      cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }

  string s = toSAM(r,reader);
  if(!gzwrite(gzoutSAM,s.c_str(),s.size()))
    {
      cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }
}

void outputM( gzFile gzout,
	      gzFile gzoutSAM,
	      const bamrecord & r,
	      const bamreader & reader)
{
  assert( ! r.flag().query_unmapped );
  assert( ! r.flag().mate_unmapped );
  string name = editRname(r.read_name());
  vector<mapping_pos> mpos = get_mapping_pos(r,reader);
  for( unsigned i=0;i<mpos.size();++i)
    {
      ostringstream out;
      out << name << '\t' 
	  << r.mapq() << '\t'
	  << mpos[i].chrom << '\t'
	  << mpos[i].start << '\t'
	  << mpos[i].stop << '\t'
	  << mpos[i].strand << '\t'
	  << mpos[i].mm << '\t'
	  << mpos[i].gap << '\n';

      if(!gzwrite(gzout,out.str().c_str(),out.str().size()))
	{
	  cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
    }
  string s = toSAM(r,reader);
  if(!gzwrite(gzoutSAM,s.c_str(),s.size()))
    {
      cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }
}

void updateBucket( readbucket & rb, bamrecord & b, 
		   gzFile csvfile, gzFile samfile,
		   const char * maptype,
		   const bamreader & reader )
{
  string n = editRname(b.read_name());
  auto i = rb.find(n);
  if(i == rb.end())
    {
      rb.insert( make_pair(n, std::move(b)) );
    }
  else //We've got our pair, so write it out
    {
      ostringstream o;
      if(i->second.refid() > reader.ref_cend()-reader.ref_cbegin())
	{
	  cerr << "Error: reference ID number : "<< i->second.refid()
	       << " is not present in the BAM file header. "
	       << " Line " << __LINE__ << " of " << __FILE__ << '\n';
	  exit(1);
	}
      if(b.refid() > reader.ref_cend()-reader.ref_cbegin())
	{
	  cerr << "Error: reference ID number : "<< b.refid()
	       << " is not present in the BAM file header. "
	       << " Line " << __LINE__ << " of " << __FILE__ << '\n';
	  exit(1);
	}
      auto REF = reader.ref_cbegin()+i->second.refid();
      o << editRname(i->second.read_name()) << '\t'
	<< REF->first << '\t'
	<< i->second.mapq() << '\t'
	<< i->second.pos() << '\t'
	<< i->second.pos() + alignment_length(i->second) - 1 << '\t'
	<< i->second.flag().qstrand << '\t'
	<< mismatches(i->second) << '\t'
	<< ngaps(i->second) << '\t'
	<< maptype << '\t';
      REF = reader.ref_cbegin()+b.refid();
      //Second read data
      o << editRname(b.read_name()) << '\t'
	<< REF->first << '\t'
	<< b.mapq() << '\t'
	<< b.pos() << '\t'
	<< b.pos() + alignment_length(b) - 1 << '\t'
	<< b.flag().qstrand << '\t'
	<< mismatches(b) << '\t'
	<< ngaps(b) << '\t'
	<< maptype << '\n';
      if(! gzwrite( csvfile,o.str().c_str(),o.str().size() ) )
	{
	  cerr << "Error: gzwrite error at line "
	       << __LINE__ << " of file "
	       << __FILE__ << '\n';
	  exit(1);
	}
      string SAM = toSAM(b,reader);
      if(!gzwrite(samfile,SAM.c_str(),SAM.size()))
	{
	  cerr << "Error: gzwrite error at line "
	       << __LINE__ << " of " << __FILE__ << '\n';
	  exit(1);
	}
      SAM = toSAM(i->second,reader);
      if(!gzwrite(samfile,SAM.c_str(),SAM.size()))
	{
	  cerr << "Error: gzwrite error at line "
	       << __LINE__ << " of " << __FILE__ << '\n';
	  exit(1);
	}
      rb.erase(i);
    }
}

string toSAM(const bamrecord & b,
	     const bamreader & reader)
{
  if(b.refid() > reader.ref_cend()-reader.ref_cbegin())
    {
      cerr << "Error: refID " << b.refid()
	   << " does not exist in BAM header. Line "
	   << __LINE__ << " of " << __FILE__ << '\n';
    }
  if(b.next_refid() > reader.ref_cend()-reader.ref_cbegin())
    {
      cerr << "Error: next_refID " << b.next_refid()
	   << " does not exist in BAM header. Line "
	   << __LINE__ << " of " << __FILE__ << '\n';
    }
  ostringstream buffer;
  samflag f(b.flag());
  bool um = f.query_unmapped,mum=f.mate_unmapped;
  buffer << b.read_name() << '\t'
	 << b.flag() << '\t'
	 << (!um ? (reader.ref_cbegin()+b.refid())->first : std::string("*"))<< '\t'
	 << b.pos()+1 << '\t'
	 << b.mapq() << '\t'
	 << b.cigar() << '\t'
	 << (!mum ? (reader.ref_cbegin()+b.next_refid())->first : std::string("*") )<< '\t'
	 << b.next_pos()+1 << '\t'
	 << ((!um && !mum) ? b.tlen() : 0)<< '\t'
	 << b.seq() << '\t';
  std::for_each(b.qual_cbeg(),b.qual_cend(),[&](const char & ch) {
      buffer << char(ch+33);
    });
  buffer << '\t'
	 << b.allaux() << '\n';
  return buffer.str();
}
