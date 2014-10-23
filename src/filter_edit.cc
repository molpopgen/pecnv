#include <map>
#include <utility>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <cassert> 
#include <boost/bind.hpp>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <sstream>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <cctype>
#include <isbinary.hpp>

using namespace std;
using namespace boost;
using namespace boost::iostreams;

void read_file (const char * file, char * outfile, int filter);

int main(int argc, char ** argv)
{
  
  int argn = 1;
  const char * file = argv[argn++];
  char * outfile = argv[argn++];//for novel
  const int filter = atoi(argv[argn++]);
  
  read_file(file, outfile, filter);
  
}

template<typename streamtype> 
void printout( streamtype & in, char * outfile, int filter)
{
  int nplus, nminus;
  string temp, pdist, mdist, mfirst, mlast, chromo, pfirst, plast, pin, min;

  filtering_ostream outstream;
  outstream.push(gzip_compressor());
  outstream.push(file_sink(outfile,ios_base::out|ios_base::binary));  
  
  getline(in, temp);
  int maxdist;
  string na = ("NA");
  maxdist = 1000;
  outstream << "chromo" << "\t" <<  "nplus" << "\t" << "nminus" << "\t" << "pfirst" << "\t" << "plast" << "\t" << "pdist" << "\t" << "pin" << "\t" << "mfirst" << "\t" << "mlast" << "\t" << "mdist" << "\t" << "min" << endl;

  while(! in.eof() )
    {

      int pdisttemp, mdisttemp, pintemp, mintemp;

      pdisttemp = -1;
      mdisttemp = -1;
      pintemp = -1;
      mintemp = -1;
	
      in >> chromo >> nplus >> nminus >> pfirst >> plast >> pdist >> pin >> mfirst >> mlast >> mdist >> min >> ws;
      // cout << chromo << "\t" << nplus << "\t" << nminus << "\t" << pfirst << "\t" << plast << "\t" << pdist << "\t" << pin << "\t" << mfirst << "\t" << mlast << "\t" << mdist << "\t" << min << endl;

      //first check to see if pdist and mdist are NA or values

      if(pdist.compare(na) != 0) {//not the same

	pdisttemp = atoi(pdist.c_str());

      }
      if(mdist.compare(na) != 0) {//not the same

	mdisttemp = atoi(mdist.c_str());

      }     
      if(pin.compare(na) != 0) {//not the same

	pintemp = atoi(pin.c_str());

      }
      if(min.compare(na) != 0) {//not the same

	mintemp = atoi(min.c_str());

      }  
       
      if( (nplus >= filter ) and (nminus >= filter)) // first filter on required number of reads
	  {
	    if((pdisttemp >= maxdist) and (mdisttemp >= maxdist) and (pintemp == 0) and (mintemp == 0)) {
	      
	      outstream << chromo << "\t" <<  nplus << "\t" << nminus << "\t" << pfirst << "\t" << plast << "\t" << pdist << "\t" << pin << "\t" << mfirst << "\t" << mlast << "\t" << mdist << "\t" << min << endl;
	    }
	  }
    }
  
}

void read_file (const char * file, char * outfile,  int filter)
{
  if( isbinary(file) )
    {
      filtering_istream in;
      in.push(gzip_decompressor());
      in.push(file_source(file,ios_base::in|ios_base::binary));
      printout(in,outfile,filter);
    }
  else
    {
      ifstream in(file);
      printout(in,outfile, filter);
    }
}
