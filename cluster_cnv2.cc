/*  
  cluster_cnv2.cc

  Faster version that only reads input files once.

  Copyright 2010 Kevin Thornton, University of California Irvine

  This code is released under the terms of the GNU Public License

  This code reads in the *_structural*.csv.gz from a single line,
  and does the following:

  For each chromosome, read pairs are clustered into the same
  CNV call if the following criteria are met:

  1. Two reads have the same lane and read id
  OR
  2. Reads from different read pairs map to the same chromosome
  within mdist base pairs of each other on the same strand
*/

#include <cstdlib>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <cassert>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <limits>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/bind.hpp>

#include <string_unsigned_lookup.hpp>
#include <isbinary.hpp>
#include <file_util.hpp>

using namespace std;
using namespace boost;
using namespace boost::iostreams;

const unsigned UMAX = numeric_limits<unsigned>::max();
std::string make_readname(const unsigned & line,
			  const unsigned & lane, const unsigned & pair);

struct linkeddata
{
  mutable unsigned a,aS,b,bS; //positions on strands -- start1,stop1,start2,stop2
  mutable string readname;
  short strand1,strand2;
  linkeddata(const unsigned & __a, 
	     const unsigned & __aS, 
	     const unsigned & __b,
	     const unsigned & __bS,
	     const string & __readname,
	     /*
	       const unsigned & line,
	       const unsigned & lane, const unsigned & pair,
	     */
	     const short & _strand1,
	     const short & _strand2) : a(__a),
				       aS(__aS),
				       b(__b),
				       bS(__bS),
				       readname( __readname ), //make_readname(line,lane,pair) ),
				       strand1(_strand1),strand2(_strand2)
  {
  }
};

struct sort_linked : public binary_function<linkeddata,linkeddata,bool>
{
  inline bool operator()(const linkeddata & lhs,const linkeddata & rhs) const
  {
    return lhs.a < rhs.a && lhs.b < rhs.b;
  }
};

struct order_clusters : public binary_function< vector<vector<linkeddata>::const_iterator>,
						vector<vector<linkeddata>::const_iterator>,
						bool >
{
  inline bool operator()( const vector<vector<linkeddata>::const_iterator> & a,
			  const vector<vector<linkeddata>::const_iterator> & b ) const
  {
    unsigned min1 = numeric_limits<unsigned>::max(),
      min2 = numeric_limits<unsigned>::max();
    for(unsigned i=0;i<a.size();++i)
      {
	min1=min(min1,a[i]->a);
      }
    for(unsigned i=0;i<b.size();++i)
      {
	min2=min(min2,b[i]->a);
      }
    return min1 < min2;
  };
};

typedef vector< vector<vector<linkeddata>::const_iterator> >  cluster_container;


cluster_container cluster_linked( const vector<linkeddata> & raw,
				  const unsigned & mdist );

void reduce_clusters( cluster_container & clusters,
		      const unsigned & mdist );

unsigned mindist(const unsigned & st1,
		 const unsigned & stp1,
		 const unsigned & st2,
		 const unsigned & stp2);

bool pair_should_cluster( vector<linkeddata>::const_iterator & pair,
			  vector<vector<linkeddata>::const_iterator> & cluster,
			  const unsigned & mdist);

bool unique_positions(const vector<linkeddata> & data,
		      const unsigned & start,
		      const unsigned & stop,
		      const unsigned & start2,
		      const unsigned & stop2);

/*
bool unique_positions_by_lane(const vector<linkeddata> & data,
			      const unsigned & start,
			      const unsigned & stop,
			      const unsigned & start2,
			      const unsigned & stop2,
			      const unsigned & lane,
			      const bool & check_lane = false);
*/

void write_clusters( filtering_ostream & o,
		     const string & chrom1,
		     const string & chrom2,
		     const cluster_container & clusters,
		     unsigned * eventid );


template<typename streamtype>
void read_data_details(map<unsigned,vector<linkeddata> > & raw_div,
		       map<unsigned,vector<linkeddata> > & raw_par,
		       map<unsigned,map<unsigned,vector<linkeddata> > > & raw_ul,
		       streamtype & lin,
		       streamtype & rin,	      
		       const unsigned & min_mqual,
		       const unsigned & max_mm,
		       const unsigned & max_gap,
		       vector<pair<string,unsigned> > * chrom_labels,
		       unsigned * chrom_index)
{
  string chrom_label,chrom_label2,pairname,pairname2;
  //unsigned line,lane,pair,read,
  unsigned mqual,chrom,strand,mm,gap,
    //line2,lane2,pair2,read2,
    mqual2,chrom2,strand2,mm2,gap2;

  int start,stop,start2,stop2;
  std::string type,type2;

  while(! lin.eof() )
    {
      /*
	lin >> line >> lane >> pair >> read 
	>> mqual >> chrom_label >> start >> stop >> strand >> mm >> gap >> type >> ws;
	rin >> line2 >> lane2 >> pair2 >> read2 
	>> mqual2 >> chrom_label2 >> start2 >> stop2 >> strand2 >> mm2 >> gap2 >> type2 >> ws;
      */
      lin >> pairname
	  >> mqual >> chrom_label >> start >> stop >> strand >> mm >> gap >> type >> ws;
      rin >> pairname2
	  >> mqual2 >> chrom_label2 >> start2 >> stop2 >> strand2 >> mm2 >> gap2 >> type2 >> ws;
      assert(pairname == pairname2);
      chrom = update_lookup(chrom_labels,chrom_index,chrom_label);
      chrom2 = update_lookup(chrom_labels,chrom_index,chrom_label2);
      /*
      assert(type==type2);
      assert(line==line2);
      assert(lane==lane2);
      assert(pair==pair2);
      */
      if( mqual >= min_mqual && mqual2 >= min_mqual &&
	  mm <= max_mm && gap <= max_gap &&
	  mm2 <= max_mm && gap2 <= max_gap )
	{
	  if(type == "DIV")
	    {
	      /*
	      if ( unique_positions_by_lane(raw_div[chrom],
						(strand==0) ? start2 : start,
						(strand==0) ? stop2 : stop,
						(strand==0) ? start : start2,
						(strand==0) ? stop : stop2,
						lane) )
		    {
	      */
	      if ( unique_positions(raw_div[chrom],
				    (strand==0) ? start2 : start,
				    (strand==0) ? stop2 : stop,
				    (strand==0) ? start : start2,
				    (strand==0) ? stop : stop2 ) )
		{
		  assert( (strand==0) ? (strand2 == 1) : (strand == 1) );
		  raw_div[chrom].push_back( linkeddata( (strand==0) ? start2 : start,
							(strand==0) ? stop2 : stop,
							(strand==0) ? start : start2,
							(strand==0) ? stop : stop2,
							pairname,1,0 ) );
		}
						    //line,lane,pair,1,0 ) );
		      /*
			}
		      */
	    }
	  else if (type == "PAR")
	    {
	      /*
	      if ( unique_positions_by_lane(raw_par[chrom],
					    (start<start2) ? start : start2,
					    (start<start2) ? stop : stop2,
					    (start<start2) ? start2 : start,
					    (start<start2) ? stop2 : stop,
					    lane) )
		{
	      */
	      if ( unique_positions(raw_par[chrom],
				    (start<start2) ? start : start2,
				    (start<start2) ? stop : stop2,
				    (start<start2) ? start2 : start,
				    (start<start2) ? stop2 : stop) )
		{
		  raw_par[chrom].push_back( linkeddata( (start<start2) ? start : start2,
							(start<start2) ? stop : stop2,
							(start<start2) ? start2 : start,
							(start<start2) ? stop2 : stop,
							pairname,//line,lane,pair,
							(start<start2) ? strand : strand2,
							(start<start2) ? strand2 : strand ));
		}
		  /*
		    }
		  */
	    }
	  else if (type == "UL")
	    {
	      assert(chrom != chrom2);
	      if( chrom > chrom2 )
		{
		  swap(chrom,chrom2);
		  swap(start,start2);
		  swap(stop,stop2);
		  swap(strand,strand2);
		}
	      /*
		if ( unique_positions_by_lane(raw_ul[chrom][chrom2],start,stop,start2,stop2,lane) )
		{
	      */
	      if ( unique_positions(raw_ul[chrom][chrom2],start,stop,start2,stop2) )
		{
		  raw_ul[chrom][chrom2].push_back( linkeddata(start,stop,start2,stop2,
							      pairname,strand,strand2) );
		}
	      //line,lane,pair,strand,strand2) );
	      /*
		}
	      */
	    }
#ifndef NDEBUG
	  else
	    {
	      abort();
	    }
#endif
	}
    }
}

void read_data( map<unsigned,vector<linkeddata> > & raw_div,
		map<unsigned,vector<linkeddata> > & raw_par,
		map<unsigned,map<unsigned,vector<linkeddata> > > & raw_ul,
		const char * left,
		const char * right,
		const unsigned & min_mqual,
		const unsigned & max_mm,
		const unsigned & max_gap,
		vector<pair<string,unsigned> > * chrom_labels,
		unsigned * chrom_index)
{
  if( isbinary(left) && isbinary(right) )
    {
      filtering_istream il;
      il.push(gzip_decompressor());
      il.push(file_source(left,ios_base::in|ios_base::binary));
      
      filtering_istream ir;
      ir.push(gzip_decompressor());
      ir.push(file_source(right,ios_base::in|ios_base::binary));

      read_data_details( raw_div,raw_par,raw_ul , il,ir,min_mqual,max_mm,max_gap,chrom_labels,chrom_index );
    }
  else
    {
      ifstream il(left),ir(right);
      read_data_details( raw_div,raw_par,raw_ul , il,ir,min_mqual,max_mm,max_gap,chrom_labels,chrom_index );
    }
}
		


int main(int argc, char ** argv)
{
  int argn=1;
  if( argc < 8 )
    {
      cerr << "usage: "
	   << argv[0]
	   << " min_quality max_mm max_gap insert_size outfile_div outfile_par outfile_ul "
	   << "structural_file1a structural_file1b ... structural_fileNa structural_fileNb\n";
      exit(0);
    }
  const unsigned min_mqual = atoi(argv[argn++]);
  const unsigned max_mm = atoi(argv[argn++]);
  const unsigned max_gap = atoi(argv[argn++]);
  const unsigned mdist = atoi(argv[argn++]);
  const char * divfile = argv[argn++];
  const char * parfile = argv[argn++];
  const char * ulfile = argv[argn++];

  //make sure output files are writable
  const string header = "id\tchrom1\tcoverage\tstrand1\tstart1\tstop1\tchrom2\tstrand2\tstart2\tstop2\treads";

  if( file_exists(divfile) )
    {
      cerr << "error: " << divfile << " already exists, and we don't want to accidentally over-write something important!\n";
      exit(10);
    }
  filtering_ostream divstream;
  divstream.push(gzip_compressor());
  divstream.push(file_sink(divfile,ios_base::out|ios_base::binary));
  divstream << header << '\n';

  if( file_exists(parfile) )
    {
      cerr << "error: " << parfile << " already exists, and we don't want to accidentally over-write something important!\n";
      exit(10);
    }
  filtering_ostream parstream;
  parstream.push(gzip_compressor());
  parstream.push(file_sink(parfile,ios_base::out|ios_base::binary));
  parstream << header << '\n';

  if( file_exists(ulfile) )
    {
      cerr << "error: " << ulfile << " already exists, and we don't want to accidentally over-write something important!\n";
      exit(10);
    }
  filtering_ostream ulstream;
  ulstream.push(gzip_compressor());
  ulstream.push(file_sink(ulfile,ios_base::out|ios_base::binary));
  ulstream << header << '\n';

  map<unsigned, vector<linkeddata> > raw_div;
  map<unsigned, vector<linkeddata> > raw_par;
  map<unsigned, map<unsigned,vector<linkeddata> > > raw_ul;
  vector<pair<string,unsigned> > chrom_labels;
  unsigned chrom_index=0;
  for(int i = argn;i<argc;i+=2)
    {
      cerr << "processing " << argv[i]
	   << " and " << argv[i+1] << '\n';
      read_data(raw_div,raw_par,raw_ul,
		argv[i],argv[i+1],min_mqual,max_mm,max_gap,&chrom_labels,&chrom_index);
    }

  unsigned eventid=0;
  cerr << "clustering div\n";
  for(map<unsigned,vector<linkeddata> >::iterator itr = raw_div.begin();
      itr != raw_div.end();++itr)
    {
      sort(itr->second.begin(),
	   itr->second.end(),sort_linked());
      cluster_container clusters = cluster_linked(itr->second,mdist);
      cluster_container clusters2(clusters);
      sort(clusters.begin(),clusters.end(),order_clusters());
      write_clusters( divstream, 
		      lookup_string(chrom_labels,itr->first), 
		      lookup_string(chrom_labels,itr->first),
		      clusters,&eventid );
    }

  cerr << "clustering par\n";
  eventid=0;
  for(map<unsigned,vector<linkeddata> >::iterator itr = raw_par.begin();
      itr != raw_par.end();++itr)
    {
      sort(itr->second.begin(),
	   itr->second.end(),sort_linked());
      cluster_container clusters = cluster_linked(itr->second,mdist);
      sort(clusters.begin(),clusters.end(),order_clusters());
      write_clusters( parstream, 
		      lookup_string(chrom_labels,itr->first),
		      lookup_string(chrom_labels,itr->first),
		      clusters,&eventid );
    }

  cerr << "clustering ul\n";
  eventid=0;
  for( map<unsigned, map<unsigned,vector<linkeddata> > >::iterator itr = raw_ul.begin() ;
       itr != raw_ul.end() ; ++itr )
    {
      for( map<unsigned,vector<linkeddata> >::iterator itr2 = itr->second.begin() ; 
	   itr2 != itr->second.end() ; ++itr2 )
	{
	  assert(itr->first < itr2->first);
	  sort(itr2->second.begin(),itr2->second.end(),sort_linked());
	  cluster_container clusters = cluster_linked(itr2->second,mdist);
	  sort(clusters.begin(),clusters.end(),order_clusters());
	  write_clusters( ulstream, 
			  lookup_string(chrom_labels,itr->first),
			  lookup_string(chrom_labels,itr2->first),
			  clusters,&eventid );
	}
    }
}

bool unique_positions(const vector<linkeddata> & data,
		      const unsigned & start,
		      const unsigned & stop,
		      const unsigned & start2,
		      const unsigned & stop2)
{
  for(unsigned i=0;i<data.size();++i)
    {
      if( start == data[i].a && 
	  start2 == data[i].b )
	{
	  return false;
	}
    }
  return true;
}

/*
bool unique_positions_by_lane(const vector<linkeddata> & data,
			      const unsigned & start,
			      const unsigned & stop,
			      const unsigned & start2,
			      const unsigned & stop2,
			      const unsigned & lane,
			      const bool & check_lane)
{
  for(unsigned i=0;i<data.size();++i)
    {
      if( start == data[i].a &&
	  //stop == data[i].aS &&
	  start2 == data[i].b )
	  //stop2 == data[i].bS )
	{
	  if(check_lane)
	    {
	      string temp(data[i].readname);
	      string::size_type colon1 = temp.find(':');
	      string::size_type colon2 = temp.find(':',colon1+1);
	      unsigned this_read_pair_lane = atoi( string(temp.begin()+colon1+1,temp.begin()+colon2).c_str() );
	      if( lane == this_read_pair_lane )
		{
		  return false;
		}
	    }
	  else
	    {
	      return false;
	    }
	}
    }
  return true;
}
*/

/*
  OLD VERSION, USED IN DGRP
bool pair_should_cluster( vector<linkeddata>::const_iterator & pair,
			  vector<vector<linkeddata>::const_iterator> & cluster,
			  const unsigned & mdist)
{
  for( unsigned i = 0 ; i < cluster.size() ; ++i )
    {
      if( pair->strand1 == cluster[i]->strand1
	  && pair->strand2 == cluster[i]->strand2 )
	{
	  if ( (max(pair->a,cluster[i]->a)-min(pair->a,cluster[i]->a) <= mdist || //a1 vs a2
		max(pair->a,cluster[i]->aS)-min(pair->a,cluster[i]->aS) <= mdist ||// a1 vs aS2
		max(pair->aS,cluster[i]->a)-min(pair->aS,cluster[i]->a) <= mdist ||// aS1 vs a2
		max(pair->aS,cluster[i]->aS)-min(pair->aS,cluster[i]->aS) <= mdist )//as1 vs aS2
	       &&
	       (max(pair->b,cluster[i]->b)-min(pair->b,cluster[i]->b) <= mdist ||//b1 vs b2
		max(pair->b,cluster[i]->bS)-min(pair->b,cluster[i]->bS) <= mdist ||//b1 vs bS2
		max(pair->bS,cluster[i]->b)-min(pair->bS,cluster[i]->b) <= mdist ||//bS vs b2
		max(pair->bS,cluster[i]->bS)-min(pair->bS,cluster[i]->bS) <= mdist ) )//bS vs bS2
	    {
	      return true;
	    }
	}
      else
	{
	  return false;
	}
    }
  return false;
}
*/

unsigned mindist(const unsigned & st1,
		 const unsigned & stp1,
		 const unsigned & st2,
		 const unsigned & stp2)
{
  unsigned a = max(st1,st2) - min(st1,st2);
  unsigned b = max(st1,stp2) - min(st1,stp2);
  unsigned c = max(stp1,stp2) - min(stp1,stp2);
  unsigned d = max(stp1,st2) - min(stp1,st2);

  return( min(a,min(min(b,c),d)) );
}



bool pair_should_cluster( vector<linkeddata>::const_iterator & pair,
			  vector<vector<linkeddata>::const_iterator> & cluster,
			  const unsigned & mdist)
{
  for( unsigned i = 0 ; i < cluster.size() ; ++i )
    {
      if( pair->strand1 == cluster[i]->strand1
	  && pair->strand2 == cluster[i]->strand2 )
	{
	  if( mindist(pair->a,pair->aS,cluster[i]->a,cluster[i]->aS) <= mdist &&
	      mindist(pair->b,pair->bS,cluster[i]->b,cluster[i]->bS) <= mdist )
	    {
	      return true;
	    }
	}
      else if(pair->strand1 == cluster[i]->strand2 &&
	      pair->strand2 == cluster[i]->strand1)
	{
	  if( mindist(pair->a,pair->aS,cluster[i]->b,cluster[i]->bS) <= mdist &&
	      mindist(pair->b,pair->bS,cluster[i]->a,cluster[i]->aS) <= mdist )
	    {
	      return true;
	    }
	}
    }
  return false;
}

cluster_container cluster_linked( const vector<linkeddata> & raw,
				  const unsigned & mdist )
{
  typedef vector<linkeddata>::const_iterator citr;
  cluster_container clusters;

  clusters.push_back( vector<citr>(1,raw.begin()) );
  for( citr i = raw.begin()+1; i < raw.end();++i )
    {
      bool clustered=false;
      for(unsigned j=0;j<clusters.size();++j)
	{
	  if( pair_should_cluster(i,clusters[j],mdist) )
	    {
	      clusters[j].push_back(i);
	      clustered=true;
	      j=clusters.size();
	    } 
	}
      if(!clustered)
	{
	  clusters.push_back( vector<citr>(1,i) );
	}
    }
  reduce_clusters(clusters,mdist);
  return clusters;
}

void reduce_clusters( cluster_container & clusters,
		      const unsigned & mdist )
{
  typedef cluster_container::iterator citr;
  citr i=clusters.end()-1,j,beg=clusters.begin();
  while( i>beg )
    {
      bool merged = false;
      for( j = i - 1 ; !merged&&j >= beg ; --j )
	{
	  //do any reads in i cluster with any data in j?
	  for(unsigned k=0 ; !merged && k < i->size() ; ++k )
	    {
	      if( pair_should_cluster(*(i->begin()+k),*j,mdist ) )
		{
		  //answer is yes, so we merge i into j, 
		  //delete i, and take care of iterator invalidation
		  copy( i->begin(),i->end(),std::back_inserter(*j) );
		  clusters.erase(i);
		  i=clusters.end()-1;
		  beg=clusters.begin();
		  merged=true;
		}
	    }
	}
      if(!merged)
	{
	  --i;
	}
    }
}

void write_clusters( filtering_ostream & o,
		     const string & chrom1,
		     const string & chrom2,
		     const cluster_container & clusters,
		     unsigned * eventid )
{
  for(unsigned i=0;i<clusters.size();++i)
    {
      //get the boundaries of each event
      unsigned min1=UMAX,max1=0,min2=UMAX,max2=0;
      string readnames;
      for(unsigned j=0;j<clusters[i].size();++j)
	{
	  min1 = min(min1,clusters[i][j]->a);
	  max1 = max(max1,clusters[i][j]->aS);
	  min2 = min(min2,clusters[i][j]->b);
	  max2 = max(max2,clusters[i][j]->bS);
	  ostringstream t;
	  t << ';' 
	    << clusters[i][j]->a << ',' 
	    << clusters[i][j]->aS << ','
	    << clusters[i][j]->strand1 << ','
	    << clusters[i][j]->b << ',' 
	    << clusters[i][j]->bS << ','
	    << clusters[i][j]->strand2;
	  if ( readnames.empty() )
	    {
	      readnames += clusters[i][j]->readname;
	    }
	  else
	    {
	      readnames += "|";
	      readnames += clusters[i][j]->readname;
	    }
	  readnames += t.str();
	}
      o << *eventid << '\t'
	<< chrom1 << '\t'
	<< clusters[i].size() << '\t'
	<< clusters[i][0]->strand1 << '\t'
	<< min1 << '\t' << max1 << '\t'
	<< chrom2 << '\t'
	<< clusters[i][0]->strand2 << '\t'
	<< min2 << '\t' << max2 << '\t'
	<< readnames << '\n';
      ++(*eventid);
    } 
}

std::string make_readname(const unsigned & line,
			  const unsigned & lane, const unsigned & pair)
{
  ostringstream o;
  o << line << ':' << lane << ':' << pair;
  return o.str();
}  

