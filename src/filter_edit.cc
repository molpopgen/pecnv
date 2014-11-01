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
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <sstream>
#include <cctype>

#include <Sequence/IOhelp.hpp>

using namespace std;

void printout( const char * infile, const char * outfile, int filter);

int main(int argc, char ** argv)
{
  
  int argn = 1;
  const char * file = argv[argn++];
  char * outfile = argv[argn++];//for novel
  const int filter = atoi(argv[argn++]);
  printout(file,outfile,filter);
}

void printout( const char * infile, const char * outfile, int filter)
{
  int nplus, nminus;
  string temp, pdist, mdist, mfirst, mlast, chromo, pfirst, plast, pin, min;

  ostringstream outstream;

  gzFile gzin = gzopen(infile,"r");
  if(gzin == NULL) {
    cerr << "Error: could not open " 
	 << infile
	 << " for reading.\n";
    exit(1);
  }

  auto line = Sequence::IOhelp::gzreadline(gzin);
  int maxdist;
  string na = ("NA");
  maxdist = 1000;
  outstream << "chromo" << "\t" <<  "nplus" << "\t" << "nminus" << "\t" << "pfirst" << "\t" << "plast" << "\t" << "pdist" << "\t" << "pin" << "\t" << "mfirst" << "\t" << "mlast" << "\t" << "mdist" << "\t" << "min" << endl;

  do
    {
      int pdisttemp, mdisttemp, pintemp, mintemp;

      pdisttemp = -1;
      mdisttemp = -1;
      pintemp = -1;
      mintemp = -1;
      line = Sequence::IOhelp::gzreadline(gzin);
      istringstream in(line.first);
      in >> chromo >> nplus >> nminus >> pfirst >> plast >> pdist >> pin >> mfirst >> mlast >> mdist >> min >> ws;

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
    }while(!gzeof(gzin));
  
  gzclose(gzin);

  //Re-use the FH for output
  gzin = gzopen(outfile,"w");
  if(gzin == NULL) {
    cerr << "Error: could not open "
	 << outfile
	 << " for writing.\n";
    exit(10);
  }
  //Make sure the buffer is adequate
  //8912 is the zlib default
  gzbuffer(gzin,max(8192u,unsigned(outstream.str().size())+1));
  if( gzprintf(gzin,"%s",outstream.str().c_str()) <= 0 )
    {
      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }
  gzclose(gzin);
}

