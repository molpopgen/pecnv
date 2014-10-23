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
//pair, lane, side, chromo, event_id

typedef pair < unsigned, unsigned > puu;
typedef pair < unsigned, puu > data;
typedef map < puu, data > event;
typedef pair < int, unsigned > piu;

struct te_boundary
//store left and right together, so that we 1/2 the number of searches that need doing
{
  unsigned id,leftmin,leftmax,rightmin,rightmax;
  te_boundary() : id(0),leftmin(0),leftmax(0),rightmin(0),rightmax(0)
  {
  }
  te_boundary( const unsigned & i,
	       const unsigned & lmin,
	       const unsigned & lmax,
	       const unsigned & rmin,
	       const unsigned & rmax ) : id(i),leftmin(lmin), leftmax(lmax),rightmin(rmin),rightmax(rmax)
  {
  }
};

typedef map<string,vector<te_boundary> > tecontainer;

void read_tes( const char * teclust_datafile, const unsigned & distance, tecontainer & tes, char * outputfile2);
void read_umu ( const  tecontainer & tes, const char * umufile, event & left_events, event & right_events);
void get_mapfile ( const  tecontainer & tes, const char * mapfile, event & left_events, event & right_events);
void print_events(event & events, char * outfile);

vector < pair <int, unsigned > > overlap( const vector<te_boundary> & vbound, const unsigned & beg, const unsigned & end )
/*
  Returns a pair of values.  The first is an integer:

  0 if nothing found
  1 if overlaps a left
  2 if overlaps a right

  The second is the id of the te insertion
*/

//need to return a vector of matches vector < pair < int, unsigned > > vector will have 1 or 2 entries
{
  vector < pair < int, unsigned > > matches;
  int check = 0;
  vector<te_boundary>::const_iterator i = vbound.begin();
  
  while(i != vbound.end()) 
  {  //overlaps both right and left
      if ( (( beg >= i->leftmin && beg <= i->leftmax ) || ( end >= i->leftmin && end <= i->leftmax )) &&
	   ( ( beg >= i->rightmin && beg <= i->rightmax ) || ( end >= i->rightmin && end <= i->rightmax ) ))
	{
	  //add to vector and set 
	  check = 1;
	  matches.push_back( make_pair(3,i->id));  

	  //cout << 3 << "\t" << i->id << endl;
	}
      //overlaps left
      else if( ( beg >= i->leftmin && beg <= i->leftmax ) ||
	  ( end >= i->leftmin && end <= i->leftmax ) )
	{
	  check = 1;
	  matches.push_back(make_pair(1,i->id));  
	  // cout << 1 << "\t" << i->id << endl;
	}
      //overlaps right
     else if( ( beg >= i->rightmin && beg <= i->rightmax ) ||
	  ( end >= i->rightmin && end <= i->rightmax ) )
	{
	  check = 1;
	  matches.push_back( make_pair(2,i->id));  
	  // cout << 2 << "\t" << i->id << endl;
	}
  
      i++;
}

  if(check == 1) {

    return matches;

  }
    
  matches.push_back( make_pair(0,0u));
  
  return matches;
  
}

int main(int argc, char ** argv) 
{
  
  //if( (argc % 3) != 2 )
  // {
  //  cerr << "usage: get_pair_ids teclust_input distance outputfile outputfile2 then (um_uniqueFILE left_mappingFILE right_mappingFILE) in sets of 3 \n";
  //  exit(10);
  // }

  int argn = 1;
  const char * teclust_datafile = argv[argn++];
  const unsigned distance = atoi(argv[argn++]);
  char * outputfile = argv[argn++];
  char * outputfile2 = argv[argn++];

  tecontainer tes;

  read_tes(teclust_datafile,distance,tes,outputfile2);

  while (argn < argc) 
    {
      const char * um_unique_datafile = argv[argn++];
      //const char * left_mapping_datafile = argv[argn++];
      //const char * right_mapping_datafile = argv[argn++];

      event left_events, right_events;      

      read_umu(tes, um_unique_datafile, left_events, right_events);
      //      get_mapfile(tes, left_mapping_datafile, left_events, right_events);
      //      get_mapfile(tes, right_mapping_datafile, left_events, right_events);//these will print to the output file 
      print_events(left_events, outputfile);
      print_events(right_events, outputfile);
  }
}

void read_tes( const char * teclust_datafile, const unsigned & distance, tecontainer & tes, char * outputfile2)
{
  //read in teclust_datafile 
  cout << "reading teclust_input" << endl;
  string temp, pdist, mdist;
  string chromo;
  unsigned nplus, nminus, pfirst, plast, pin, mfirst, mlast, min;
  unsigned counter=0;
  
  ofstream outfile2;
  outfile2.open(outputfile2, ofstream::app);

  if( isbinary(teclust_datafile) )
    {
      filtering_istream in;
      in.push(gzip_decompressor());
      in.push(file_source(teclust_datafile,ios_base::in|ios_base::binary));
      getline(in, temp);
      while(!in.eof())
	{
	  in >> chromo >> nplus >> nminus >> pfirst >> plast >> pdist >> pin >> mfirst >> mlast >> mdist >> min >> ws;
	  outfile2 << chromo << "\t" << counter << "\t" << plast << "\t" <<  mfirst << endl;
	  tes[chromo].push_back( te_boundary(counter++,pfirst-distance,plast,mfirst,mlast+distance ));

	} 
	
    }
  else
    {
      ifstream file(teclust_datafile);
      getline(file, temp);
      while(!file.eof())
	{
	  file >> chromo >> nplus >> nminus >> pfirst >> plast >> pdist >> pin >> mfirst >> mlast >> mdist >> min >> ws;
	  outfile2 << chromo << "\t" << counter << "\t" << plast << "\t" << mfirst << endl;
	  tes[chromo].push_back( te_boundary(counter++, pfirst-distance, plast, mfirst, mlast+distance )) ;

	}
    }
  outfile2.close();
}

template<typename streamtype> 
void read_umu_details( streamtype & in, 
		       const  tecontainer & tes, event & left_events, event & right_events)
{
  unsigned line, lane, pair, read, mapping, beg, end, strand, mismatch, gap;
  string chromo;
  cout << "reading umu" << endl;
  while(! in.eof() )
    {
      in >> line >> lane >> pair >> read >> mapping >> chromo >> beg >> end >> strand >> mismatch >> gap >> ws;
      tecontainer::const_iterator itr = tes.find(chromo);
      if( itr != tes.end() )	
	{
	  //does read overlap a left or right?
	  vector < piu > rv = overlap( itr->second, beg,end );
	  vector < piu >::const_iterator j = rv.begin();
	  while (j != rv.end())
	    {
	      if ( (j->first == 1 && strand == 0)  )
		{ //pair, lane, side chrom, event
		  left_events.insert(make_pair(make_pair(lane, pair), make_pair(0, make_pair( itr->first, j->second))));
		}
	      if ((j->first == 2 && strand == 1) )
		{
		  right_events.insert(make_pair(make_pair(lane, pair), make_pair(1, make_pair( itr->first, j->second))));
		}
	      if (j->first == 3)
		{
		  left_events.insert(make_pair(make_pair(lane, pair), make_pair(0, make_pair( itr->first, j->second))));
		  right_events.insert(make_pair(make_pair(lane, pair), make_pair(1, make_pair( itr->first, j->second))));
		}

	      j++;
	    }
	}
    }
}

template<typename streamtype> 
void read_map_details( streamtype & in, 
		       const  tecontainer & tes,
		       event & left_events, event & right_events)
{
  unsigned line, lane, pair, read, mapping, beg, end, strand, mismatch, gap;
  string chromo;
  string type;
  cout << "reading mapping" << endl;
  while(! in.eof() )
    {
      in >> line >> lane >> pair >> read >> mapping >> chromo >> beg >> end >> strand >> mismatch >> gap >> type >> ws;
      tecontainer::const_iterator itr = tes.find(chromo);
      if( itr != tes.end() )
	{
	  //does read overlap a left or right?
	  std::vector < piu > rv = overlap( itr->second, beg,end );
	  vector < piu >::const_iterator j = rv.begin();
	  while ( j != rv.end()) 
	    {
	      if (j->first == 1 && (type == "UU") )
		{
		  left_events.insert(make_pair(make_pair(lane, pair), make_pair(0, make_pair( itr->first, j->second))));
		}
	      if (j->first == 1 && (type == "U") && strand == 0) 
		{
		  left_events.insert(make_pair(make_pair(lane, pair), make_pair(0, make_pair( itr->first, j->second))));
		}
	      
	      if (j->first == 2 && (type == "UU") ) 
		{
		  right_events.insert(make_pair(make_pair(lane, pair), make_pair(1, make_pair( itr->first, j->second))));
		}
	      if (j->first == 2 && (type == "U") && strand == 1 )
		{
		  right_events.insert(make_pair(make_pair(lane, pair), make_pair(1, make_pair( itr->first, j->second))));
		}
	      if (j->first == 3)
		{
		  left_events.insert(make_pair(make_pair(lane, pair), make_pair(0, make_pair( itr->first, j->second))));
		  right_events.insert(make_pair(make_pair(lane, pair), make_pair(1, make_pair( itr->first, j->second))));
		}
	      j++;
	    }
	}
    }
}

void read_umu ( const  tecontainer & tes, const char * umufile, event & left_events, event & right_events)
{
 if( isbinary(umufile) )
    {
      filtering_istream in;
      in.push(gzip_decompressor());
      in.push(file_source(umufile,ios_base::in|ios_base::binary));
      read_umu_details(in,tes,left_events,right_events);
    }
  else
    {
      ifstream in(umufile);
      read_umu_details(in,tes,left_events,right_events);
    }
}

void get_mapfile ( const  tecontainer & tes, const char * mapfile, event & left_events, event & right_events)
{
 if( isbinary(mapfile) )
    {
      filtering_istream in;
      in.push(gzip_decompressor());
      in.push(file_source(mapfile,ios_base::in|ios_base::binary));
      read_map_details(in,tes,left_events,right_events);
    }
  else
    {
      ifstream in(mapfile);
      read_map_details(in,tes,left_events,right_events);
    }
}

void print_events( event & events, char * outputfile)

{
  cout << "printing" << endl;
  ofstream outfile;
  
  outfile.open(outputfile, ofstream::app);

  for( event::const_iterator iter = events.begin(); iter != events.end(); ++iter) {
    //order out is chromo, event, lane, pair, side

    outfile << iter->second.second.first << "\t" << iter->second.second.second << "\t" << iter->first.first << "\t" << iter->first.second << "\t" << iter->second.first << endl;
  }
  outfile.close();
} 
