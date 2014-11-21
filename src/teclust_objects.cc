#include <teclust_objects.hpp>
#include <limits>
#include <iostream>
using namespace std;

params::params() : reference_datafile(string()),
		   outfile(string()),
		   bamfile(string()),
		   readfile(string()),
		   umufile(string()),
		   ummfile(string()),
		   phrapdir(string()),
		   INSERTSIZE(numeric_limits<unsigned>::max()),
		   MDIST(numeric_limits<unsigned>::max()),
		   MINREADS(numeric_limits<unsigned>::max()),
		   CLOSEST(-1),
		   novelOnly(true),
		   greedy(true)
{
}

int32_t teinfo::start() const { return this->first; }
int32_t teinfo::stop() const { return this->second; }

teinfo::teinfo( int32_t __s, int32_t __st ) : pair<int32_t,int32_t>(__s,__st)
{
}

cluster::cluster() : positions(make_pair(numeric_limits<unsigned>::max(),numeric_limits<unsigned>::max())),nreads(0)
  {
  }

cluster::cluster(const int32_t & pos1,
		 const int32_t & pos2,
		 const unsigned & nr) : positions(make_pair(pos1,pos2)),nreads(nr)
{
}
