#include <string_unsigned_lookup.hpp>
#include <algorithm>
#include <functional>
using namespace std;

struct find_string : public binary_function< pair<string,unsigned>,
					     string, bool >
{
  inline bool operator()(const pair<string,unsigned> & psu,
		       const string & s ) const
  {
    return psu.first == s;
  }
};

struct find_unsigned : public binary_function< pair<string,unsigned>,
					       unsigned, bool >
{
  inline bool operator()(const pair<string,unsigned> & psu,
		       const unsigned & u ) const
  {
    return psu.second == u;
  }
};

unsigned update_lookup( vector<pair<string,unsigned> > * vpsu,
			unsigned * index,
			const string & s )
{
  vector<pair<string,unsigned> >::iterator i = find_if(vpsu->begin(),vpsu->end(),bind2nd(find_string(),s));

if ( i == vpsu->end() )
  {
    vpsu->push_back(make_pair(s,*index));
    ++(*index);
    return ((*index)-1);
  }
 
 return i->second;
}

unsigned lookup_unsigned( const std::vector<std::pair<std::string,unsigned> > & vpsu,
			  const std::string & s )
{
  vector<pair<string,unsigned> >::const_iterator i = find_if(vpsu.begin(),vpsu.end(),bind2nd(find_string(),s));
  return i->second;
}

string lookup_string( const vector<pair<string,unsigned> > & vpsu,
	       const unsigned & u )
{
  vector<pair<string,unsigned> >::const_iterator i = find_if(vpsu.begin(),vpsu.end(),bind2nd(find_unsigned(),u));
return i->first;
}
