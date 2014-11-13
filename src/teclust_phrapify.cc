#include <teclust_phrapify.hpp>
#include <sstream>
#include <iostream>
#include <vector>

#include <Sequence/bamreader.hpp>
#include <Sequence/samfunctions.hpp>

using namespace std;
using namespace Sequence;

struct clusteredEvent
{
  string chromo;
  int nplus,nminus,pfirst,plast,pdist,pin,mfirst,mlast,mdist,min;
  std::istream & read(istream &);
};

istream & clusteredEvent::read(istream & in)
{
  in >> chromo >> nplus >> nminus >> pfirst >> plast
      >> pdist >> pin >> mfirst >> mlast 
      >> mdist >> min;
  return in;
}

istream & operator>>(istream & in, clusteredEvent & ce)
{
  return ce.read(in);
}

vector<clusteredEvent> parseClusters(const string & clusters)
{
  vector<clusteredEvent> cEs;
  istringstream in(clusters);
  string temp;
  getline(in,temp);//get rid of header
  while(!in.eof()) //read records
    {
      clusteredEvent e;
      in >> e >> ws;
      cEs.push_back(e);
    }
  return cEs;
}

void phrapify( const params & pars,
	       const string & clusters )
{
  if(pars.phrapdir.empty()) return;
  if(pars.bamfile.empty()) return;

  vector<clusteredEvent> cEs = parseClusters(clusters);
 
}
