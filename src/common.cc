#include <common.hpp>

using namespace std;

string editRname( const string & readname )
{
  auto pound = readname.find('#');
  if(pound == string::npos) return readname;
  return string(readname.begin(),readname.begin()+pound);
}

