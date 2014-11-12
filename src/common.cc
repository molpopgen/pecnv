#include <common.hpp>
//#include <algorithm>

using namespace std;

string editRname( const string & readname )
{
  auto pound = readname.find('#');
  if(pound == string::npos) return readname;
  return string(readname.begin(),readname.begin()+pound);
}

// unsigned alen(const vector< pair<char,unsigned> > & cigar_data)
// {
//   unsigned sum=0;
//   for(unsigned i=0;i<cigar_data.size();++i)
//     {
//       if ( cigar_data[i].first == 'M' ||
// 	   cigar_data[i].first == 'I' ||
// 	   cigar_data[i].first == 'D' ||
// 	   cigar_data[i].first == 'N' )
// 	{
// 	  sum += cigar_data[i].second;
// 	}
//     }
//   return sum;
// }

// unsigned mm(const unsigned & nm,  const vector< pair<char, unsigned> > & cigar_data)
// {
//   return nm - ngaps(cigar_data);
// }

// unsigned ngaps(const vector< pair<char,unsigned> > & cigar_data)
// {
//   unsigned sum=0;
//   for(unsigned i=0;i<cigar_data.size();++i)
//     {
//       if ( cigar_data[i].first == 'I' ||
// 	   cigar_data[i].first == 'D' )
// 	{
// 	  sum += cigar_data[i].second;
// 	}
//     }
//   return sum;
// }

// vector<pair<char,unsigned> > parse_cigar(const string & cigar)
// {
//   vector<pair<char,unsigned> > cigar_data;
  
//   string::const_iterator pibeg = find_if( cigar.begin(),cigar.end(),::isdigit);
//   string::const_iterator piend = find_if( pibeg+1,cigar.end(),::isalpha );
//   char * endptr;
//   while( pibeg != cigar.end() )
//     {
//       cigar_data.push_back( make_pair( *piend, 
// 				       strtoul( string(pibeg,piend).c_str(),&endptr,10 ) ) );
//       pibeg = find_if( piend+1, cigar.end(), ::isdigit );
//       piend = find_if( pibeg+1,cigar.end(),::isalpha);
//     }
//   return cigar_data;
// }
