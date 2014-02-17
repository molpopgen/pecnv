#include <bwa_util.hpp>
#include <cstdlib>

using namespace std;
using namespace Sequence;

void getLanePair(unsigned * lane,
		 unsigned * rpair,
		 const string & name)
/*
  format is:
  line:lane:read_pair:dir
*/
{
  string::size_type colon1 = name.find(":");
  string::size_type colon2 = name.find(":",colon1+1);
  string::size_type colon3 = name.find(":",colon2+1);
  *lane = atoi( string( name.begin()+colon1+1, name.begin()+colon2 ).c_str() );
  *rpair = atoi( string( name.begin()+colon2+1, name.begin()+colon3 ).c_str() );
}

bool hasXT(const samrecord & r,const string & value)
{
  for( samrecord::tag_iterator i = r.tag_begin();
       i!=r.tag_end();++i)
    {
      if(i->tag() == "XT")
	{
	  if(i->value() == value)
	    {
	      return true;
	    }
	}
    }
  return false;
}
