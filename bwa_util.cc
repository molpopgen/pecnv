#include <bwa_util.hpp>
#include <cstdlib>

using namespace std;
using namespace Sequence;

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
