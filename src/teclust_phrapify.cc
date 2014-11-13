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

vector<clusteredEvent> parseClusters(const string & clusters,const params & pars)
//Replaces functionality of filter_edit
{
  vector<clusteredEvent> cEs;
  istringstream in(clusters);
  string temp;
  getline(in,temp);//get rid of header
  while(!in.eof()) //read records
    {
      clusteredEvent e;
      in >> e >> ws;
      if(e.nplus >= pars.MINREADS && e.nminus >= pars.MINREADS) //filter on read number
	{
	  if( e.pdist >= pars.CLOSEST && e.mdist >= pars.CLOSEST && 
	      (!pars.novelOnly || (pars.novelOnly && e.pin == 0 && e.min == 0)) )
	    {
	      cEs.push_back(e);
	    }
	}
    }
  return cEs;
}

void phrapify( const params & pars,
	       const string & clusters )
{
  if(pars.phrapdir.empty()) return;
  if(pars.bamfile.empty()) return;

  vector<clusteredEvent> cEs = parseClusters(clusters,pars);
 
}
