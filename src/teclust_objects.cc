#include <teclust_objects.hpp>
#include <limits>
using namespace std;

params::params() : reference_datafile(string()),
		   outfile(string()),
		   bamfile(string()),
		   readfile(string()),
		   umufile(string()),
		   ummfile(string()),
		   phrapdir(string()),
		   INSERTSIZE(numeric_limits<unsigned>::max()),
		   MDIST(numeric_limits<unsigned>::max())
{
}

unsigned teinfo::start() const { return this->first; }
unsigned teinfo::stop() const { return this->second; }

teinfo::teinfo( unsigned __s, unsigned __st ) : pair<unsigned,unsigned>(move(__s),move(__st))
{
}

cluster::cluster() : positions(make_pair(numeric_limits<unsigned>::max(),numeric_limits<unsigned>::max())),nreads(0)
  {
  }

cluster::cluster(const unsigned & pos1,
		 const unsigned & pos2,
		 const unsigned & nr) : positions(make_pair(pos1,pos2)),nreads(nr)
{
}
