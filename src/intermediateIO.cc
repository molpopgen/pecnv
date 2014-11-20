#include <intermediateIO.hpp>
#include <Sequence/samfunctions.hpp>
#include <cstring>

using namespace std;
using namespace Sequence;

int gzwriteCstr(gzFile of, const string & s)
{
  int rv1 = gzputs(of,s.c_str()); //this does NOT write the \0
  if( rv1 <= 0 ) return rv1;
  int c = int('\0');
  int rv2 = gzputc(of,c);
  if(rv2==-1) return rv2;
  return rv1+rv2;
}

std::pair<std::string,int> gzreadCstr( gzFile in)
{
  string s;
  int read = 0;
  int ch;
  do
    {
      ch = gzgetc(in);
      if( ch < 0 ) return make_pair(s,ch);
      if(char(ch)=='\0')break;
      s += char(ch);
      read += 1;
    }
  while(char(ch) != '\0');
  return make_pair(s,read);
}

alnInfo::alnInfo( const bamrecord & b ) : start( b.pos() ),
					  stop ( b.pos() + alignment_length(b) - 1),
					  mapq( int8_t(b.mapq()) ),
					  strand( int8_t(b.flag().qstrand) ),
					  mm ( uint16_t(mismatches(b)) ),
					  ngap( uint16_t(ngaps(b)) )
{
}

alnInfo::alnInfo( gzFile in) : start(-1),
			     stop(-1),
			     mapq(-1),
			     strand(-1),
			     mm( -1 ),
			     ngap( -1 )
{
  int rv = gzread(in,&start,sizeof(int32_t));
  if( rv <= 0 ) return;
  rv = gzread(in,&stop,sizeof(int32_t));
  if( rv <= 0 ) return;
  rv = gzread(in,&mapq,sizeof(int8_t));
  if( rv <= 0 ) return;
  rv = gzread(in,&strand,sizeof(int8_t));
  if( rv <= 0 ) return;
  rv = gzread(in,&mm,sizeof(int16_t));
  if( rv <= 0 ) return;
  rv = gzread(in,&ngap,sizeof(int16_t));
  if( rv <= 0 ) return;
}

alnInfo::alnInfo( const int32_t & __start,
		  const int32_t & __stop,
		  const int32_t & __mapq,
		  const int32_t & __strand,
		  const uint32_t & __mm,
		  const uint32_t & __gap) : start(__start),
					    stop(__stop),
					    mapq(int8_t(__mapq)),
					    strand(int8_t(__strand)),
					    mm(int16_t(__mm)),
					    ngap(int16_t(__gap))
{
}


int alnInfo::write( gzFile out )
{
  int rv = 0;
  int temp = gzwrite(out,&start,sizeof(int32_t));
  if( temp <= 0 ) return temp;
  rv += temp;
  temp = gzwrite(out,&stop,sizeof(int32_t));
  if( temp <= 0 ) return temp;
  rv += temp;
  temp = gzwrite(out,&mapq,sizeof(int8_t));
  if( temp <= 0 ) return temp;
  rv += temp;
  temp = gzwrite(out,&strand,sizeof(int8_t));
  if( temp <= 0 ) return temp;
  rv += temp;
  temp = gzwrite(out,&mm,sizeof(uint16_t));
  if( temp <= 0 ) return temp;
  rv += temp;
  temp = gzwrite(out,&ngap,sizeof(uint16_t));
  if( temp <= 0 ) return temp;
  rv += temp;
  return rv;
}
