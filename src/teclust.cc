/*
  Copyright 2010 Kevin Thornton, University of California, Irvine
  
  This code is released under the terms of the GNU Public Licesne

  Input looks like:
  2074202 6       0
  18633175        8       0
  561972  13      0
  21416313        10      1
  4720393 13      0

  Which is map_position chrom strand, corresponding to the uniquely-mapping reads in 
  unique/multi pairs, where the multis map to known TEs in the reference

  It is assumed that these have already been filtered on things like mapping quality,
  for example by using umm_te_finder

  It clusters the above into "events" suggesting TEs.
*/

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <functional>
#include <algorithm>
#include <limits>
#include <cassert>
#include <isbinary.hpp>
#include <string_unsigned_lookup.hpp>
#include <boost/bind.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

using namespace std;
using namespace boost::iostreams;

// DEFINITIONS OF DATA TYPES
typedef pair<unsigned,unsigned> puu;

const unsigned UMAX = std::numeric_limits<unsigned>::max();

struct cluster
/*
  A cluster is a group of unique reads that suggest the presence of a TE
*/
{
  puu positions;
  unsigned nreads;
  cluster() : positions(puu(UMAX,UMAX)),nreads(0)
  {
  }
  cluster(const unsigned & pos1,
	  const unsigned & pos2,
	  const unsigned & nr) : positions(make_pair(pos1,pos2)),nreads(nr)
  {
  }
};

// THE BELOW ARE STRUCTURES THAT ARE USED EITHER WITH THE "SORT" OR "FIND_IF" ALGORITHMS


template<typename T>
struct pair_less : public binary_function< pair<T,T>, pair<T,T>, bool>
{
  /*
    used to sort data--for an pair, return true if the first member of
    pair1 <= the first member of pair 2
  */
  inline bool operator()(  const pair<T,T> & lhs, const  pair<T,T> & rhs ) const
  {
    return lhs.first <= rhs.first;// && lhs.first <= lhs.second;
  }
};

struct puu_less : public binary_function<puu,puu,bool>
/*
  return true if the first element of the first puu <= the first element of the second puu
*/
{
  inline bool operator()( const puu & lhs, const puu & rhs) const
  {
    return lhs.first < rhs.first;//&& lhs.second <= rhs.second;
  }
};

struct within : public binary_function< pair<unsigned,unsigned>,unsigned, bool>
/*
  returns true if the unsigned integer is within the range specified by the pair
*/
{
  inline bool operator()(const pair<unsigned,unsigned> & lhs, const unsigned & rhs) const
  {
    return ( rhs >= lhs.first && rhs <= lhs.second );
  }
};

struct closest_plus : public binary_function< pair<unsigned,unsigned>,unsigned, bool>
/*
  used to find the closest TE in the reference 3' of a certain position
*/
{
  inline bool operator()(const pair<unsigned,unsigned> & lhs, const unsigned & rhs) const
  {
    return ( lhs.first >= rhs || ( rhs >= lhs.first && rhs <= lhs.second ) );
  }
};

struct closest_minus : public binary_function< pair<unsigned,unsigned>,unsigned, bool>
/*
  used to find the closest TE in the reference 5' of a certain position
*/
{
  inline bool operator()(const pair<unsigned,unsigned> & lhs, const unsigned & rhs) const
  {
    return ( lhs.second <= rhs || ( rhs >= lhs.first && rhs <= lhs.second ) );
  }
};

struct close_enough_minus : public binary_function< cluster,cluster,bool >
{
  inline bool operator()(const cluster & minus, const cluster & plus,
			 const unsigned & MDIST) const
  {
    if( minus.positions.first < plus.positions.second ) return false;
    if( minus.positions.first - plus.positions.second <= MDIST) return true;
    return false;
  }
};

/*
  struct sort_clusters : public binary_function< pair<cluster,cluster>,
  pair<cluster,cluster>,
  bool >
  {
  inline bool operator()(const pair<cluster,cluster> & lhs,
  const pair<cluster,cluster> & rhs) const
  {
  if( lhs.first.positions.first == UMAX ||
  rhs.first.positions.first == UMAX )
  return lhs.second.positions.first <= rhs.second.positions.first;

  if( lhs.second.positions.first == UMAX ||
  rhs.second.positions.first == UMAX )
  return lhs.first.positions.first <= rhs.first.positions.first;

       
  return (lhs.first.positions.first <= rhs.first.positions.first &&
  lhs.second.positions.first <= rhs.second.positions.first);
  }
  };
*/
//DEFINITION OF FUNCTIONS

void get_ref_te(  map< unsigned, vector< pair<unsigned,unsigned> > > & reference_te,
		  const char * reference_datafile,
		  const vector<pair<string,unsigned> > & chrom_labels );
void cluster_data( vector<pair<cluster,cluster> > & clusters,
		   const vector<puu> & raw_data, 
		   const unsigned & INSERTSIZE, const unsigned & MDIST );
void reduce_ends( vector<cluster> & clusters,
		  const unsigned & INSERTSIZE );

//void output_results( const vector<cluster> & clusters, const string & chrom_label , const vector< pair<unsigned,unsigned> > & ref_te_chromo);
void output_results( filtering_ostream & out,
		     //void output_results( ofstream & out,
		     const vector<pair<cluster,cluster> > & clusters, 
		     const string & chrom_label , const vector< pair<unsigned,unsigned> > & ref_te_chromo);
		     

template<typename streamtype>
void read_raw_data(streamtype & in, 
		   map<unsigned,vector<puu> > & raw_data,
		   vector<pair<string,unsigned> > * chrom_labels )
/*
  input data are read from STDIN, and look like:
  position chrom strand
  
*/
{
  unsigned pos,chrom,strand,dummy=0;
  string chrom_label;
  map<unsigned,vector<puu> >::iterator itr;
  unsigned nread=0;
  while(!in.eof())
    {
      in >> pos >> chrom_label >> strand  >> ws;
      ++nread;
      chrom = update_lookup(chrom_labels,&dummy,chrom_label);
      itr = raw_data.find(chrom);

      if( itr == raw_data.end() ) //new chromosome
	{
	  raw_data.insert( make_pair(chrom,
				     vector<puu>(1,puu(pos,strand))) );
	}
      else //add to existing data for chromosome
	{
	  itr->second.push_back(puu(pos,strand));
	}
    }
  cerr << nread << " lines processed\n";
  //sort data per chromosome by position
}

int main( int argc, char ** argv )
{
  if( argc < 5 )
    {
      cerr << "usage: teclust reference_datafile INSERTSIZE MDIST outfile.gz [input_files]\n";
      exit(10);
    }
  int argn = 1;
  const char * reference_datafile = argv[argn++];
  const unsigned INSERTSIZE = atoi(argv[argn++]);
  const unsigned MDIST = atoi(argv[argn++]);
  const char * outfile = argv[argn++];
  /*
    << "common\t"
    << "ttl\n";
  */
  /*
    a map is like a hash in perl
    The "key" is an unsigned integer, representing the chromosome
    The "value" is a vector of pairs of unsigned integers,
    where each pair is the start and stop of a TE in the 
    reference on the chromosome
  */
  //read the raw data in from STDIN
  map< unsigned, vector<puu> > raw_data;
  vector<pair<string,unsigned> > chrom_labels;
  for(int i = argn;i<argc;++i)
    {
      cerr << "processing " << argv[i] << '\n';
      if( isbinary(argv[i]) )
	{
	  filtering_istream in;
	  in.push(gzip_decompressor());
	  in.push(file_source(argv[i],ios_base::in|ios_base::binary));
	  if(!in)
	    {
	      cerr << "cannot read " << argv[i] << '\n';
	      exit(10);
	    }
	  read_raw_data(in,raw_data, &chrom_labels );
	}
      else
	{
	  ifstream in(argv[i]);
	  if(!in)
	    {
	      cerr << "cannot read " << argv[i] << '\n';
	      exit(10);
	    }
	  read_raw_data(in,raw_data, &chrom_labels );
	}
      for( map<unsigned,vector<puu> >::iterator itr = raw_data.begin() ; 
	   itr != raw_data.end() ; ++itr )
	{
	  std::sort(itr->second.begin(),itr->second.end(),puu_less());
	}
    }

  map< unsigned, vector< pair<unsigned,unsigned> > > reference_te;
  get_ref_te(reference_te, reference_datafile,chrom_labels);
  
  filtering_ostream out;
  out.push(gzip_compressor());
  out.push(file_sink(outfile,ios_base::out|ios_base::binary));
  
  //ofstream out(outfile);
  out << "chromo\t"
      << "nplus\t"
      << "nminus\t"
      << "pfirst\t"
      << "plast\t"
      << "pdist\t"
      << "pin\t"
      << "mfirst\t"
      << "mlast\t"
      << "mdist\t"
      << "min\n";
  //go through each chromosome, cluster results for each chromosome, and print results to screen`
  map< unsigned, vector<puu> >::iterator itr;
  for( itr = raw_data.begin() ; itr != raw_data.end() ; ++itr )
    {
      vector<pair<cluster,cluster> > clusters;
      cluster_data(clusters,itr->second,INSERTSIZE,MDIST);
      output_results(out,clusters,
		     lookup_string(chrom_labels,itr->first),
		     reference_te[itr->first]);
    }
  out.pop();
  out.pop();
  exit(0);
}

void get_ref_te(  map< unsigned, vector< pair<unsigned,unsigned> > > & reference_te,
		  const char * reference_datafile,
		  const vector<pair<string,unsigned> > & chrom_labels) 
/*
  This function reads in the start/stop positions of every TE in the reference genome
  The input file is in the following format:
  chrom start stop
  
  for each TE, and chrom is an unsigned integer representing each chromosome
*/
{
  string chrom_label;
  unsigned chrom,i,j;
  ifstream in(reference_datafile);
  if(!in)
    {
      cerr << "Could not open " << reference_datafile << " for reading\n";
      exit(10);
    }

  //Read in TE positions
  map< unsigned, vector< pair<unsigned,unsigned> > >::iterator itr;
  while(!in.eof())
    {
      in >> chrom_label >> i >> j >> ws;
      chrom = lookup_unsigned(chrom_labels,chrom_label);
      itr = reference_te.find(chrom);
      if(  itr == reference_te.end() ) //if chromosome is not present as a key
	{
	  /*
	    add this key/value -- the key is chrom, and the value is a 
	    vector<pair<unsigned,unsigned> > with 1 value, which is a pair with values i and j
	  */
	  reference_te.insert( make_pair(chrom, vector<pair<unsigned,unsigned> >(1,make_pair(i,j))) );
	}
      else //chromosome exists, so just add i and j to the vector of TE positions
	{
	  itr->second.push_back(make_pair(i,j));
	}
    }
  in.close();

  //sort TE positions by start position for each chromosome
  for( itr = reference_te.begin() ; itr != reference_te.end() ; ++itr )
    {
      sort( itr->second.begin(), itr->second.end(),puu_less() );
      //pair_less<unsigned>() );
    }
}

void cluster_data( vector<pair<cluster,cluster> > & clusters,
		   const vector<puu> & raw_data, 
		   const unsigned & INSERTSIZE, const unsigned & MDIST )
{
  vector<cluster> plus,minus;
  for( unsigned i=0;i<raw_data.size();++i )
    {
      if (raw_data[i].second == 0)
	{
	  if(plus.empty())
	    {
	      plus.push_back( cluster(raw_data[i].first,raw_data[i].first,1) );
	    }
	  else
	    {
	      bool clustered=false;
	      for(unsigned j=0;!clustered&&j<plus.size();++j)
		{
		  if( max(raw_data[i].first,plus[j].positions.second)-
		      min(raw_data[i].first,plus[j].positions.second)<= INSERTSIZE )
		    {
		      assert( raw_data[i].first >= plus[j].positions.second );
		      plus[j].positions.second = raw_data[i].first;
		      plus[j].nreads++;
		      clustered=true;
		    }

		}
	      if(!clustered)
		{
		  plus.push_back( cluster(raw_data[i].first,raw_data[i].first,1) );
		}
	      //because positions are sorted in ascending order, we only need to check for overlap
	      //with previous entry
	      /*
		if( raw_data[i].first - plus[plus.size()-1].positions.second <= INSERTSIZE )
		{
		plus[plus.size()-1].positions.second = raw_data[i].first;
		plus[plus.size()-1].nreads++;
		}
		else
		{
		plus.push_back( cluster(raw_data[i].first,raw_data[i].first,1) );
		}
	      */
	    }
	}
      else //minus strand works same as plus strand
	{
	  if( minus.empty() )
	    {
	      minus.push_back( cluster(raw_data[i].first,raw_data[i].first,1) );
	    }
	  else
	    {
	      bool clustered=false;
	      for(unsigned j=0;!clustered&&j<minus.size();++j)
		{
		  if( max(raw_data[i].first,minus[j].positions.second)-
		      min(raw_data[i].first,minus[j].positions.second) <= INSERTSIZE )
		    //if( raw_data[i].first - minus[minus.size()-1].positions.first <= INSERTSIZE )
		    {
		      assert( raw_data[i].first >= minus[j].positions.second );
		      minus[j].positions.second = raw_data[i].first;
		      minus[j].nreads++;
		      clustered=true;
		    }
		  
		}
	      if(!clustered)
		{
		  minus.push_back( cluster(raw_data[i].first,raw_data[i].first,1) );
		}
	      /*
		if( raw_data[i].first - minus[minus.size()-1].positions.second <= INSERTSIZE )
		//if( raw_data[i].first - minus[minus.size()-1].positions.first <= INSERTSIZE )
		{
		minus[minus.size()-1].positions.second = raw_data[i].first;
		minus[minus.size()-1].nreads++;
		}
		else
		{
		minus.push_back( cluster(raw_data[i].first,raw_data[i].first,1) );
		}
	      */
	    }
	}
    }
  
  reduce_ends( plus, INSERTSIZE );
  reduce_ends( minus, INSERTSIZE );
  //cerr << "here\n";
  //now, we have to match up plus and minus based on MDIST
  vector<short> matched(plus.size(),0);
  for(unsigned i=0;i<plus.size();++i)
    {
      if(!matched[i])
	{
	  vector<cluster>::iterator j = find_if(minus.begin(),
						minus.end(),
						boost::bind(close_enough_minus(),_1,plus[i],MDIST));
	  if( j != minus.end() )
	    {
	      //is there a better match in plus for this minus?
	      unsigned dist = j->positions.first-plus[i].positions.second;
	      unsigned winner = i;
	      for(unsigned k=i+1;k<plus.size();++k)
		{
		  if ( close_enough_minus()(*j,plus[k],MDIST) )
		    {
		      if( j->positions.first - plus[k].positions.second < dist )
			{
			  dist = j->positions.first - plus[k].positions.second;
			  winner=k;
			}
		    }
		}
	      clusters.push_back( make_pair(plus[winner],*j) );
	      minus.erase(j);
	      matched[winner]=1;
	      //need to take care of i, too, if i no longer matches
	      if(winner!=i)
		{
		  matched[i]=1;
		  clusters.push_back( make_pair(plus[i],cluster()) );
		}
	    }
	  else
	    {
	      matched[i]=1;
	      clusters.push_back( make_pair(plus[i],cluster()) );
	    }
	}
    }
  //cerr << "here2 " << minus.size() << '\n';
  //now, add in the minuses
  for( unsigned i = 0 ; i < minus.size() ; ++i  )
    {
      bool pushed = false;
      for(unsigned j = 0 ; j < clusters.size() ; ++j )
	{
	  if( clusters[j].first.positions.first != UMAX )
	    {
	      if( minus[i].positions.first < clusters[j].first.positions.first )
		{
		  //add it
		  pushed=true;
		  clusters.insert(clusters.begin()+j,make_pair(cluster(),minus[i]));
		  j=clusters.size();
		}
	      else if ( clusters[j].second.positions.first != UMAX )
		{
		  if( minus[i].positions.first <  clusters[j].second.positions.first )
		    {
		      pushed=true;
		      clusters.insert(clusters.begin()+j,make_pair(cluster(),minus[i]));
		      j=clusters.size();
		    }
		}
	    }
	  else if( clusters[j].second.positions.first != UMAX )
	    {
	      assert( clusters[j].second.positions.first != UMAX );
	      if( minus[i].positions.first < clusters[j].second.positions.first )
		{
		  pushed=true;
		  clusters.insert(clusters.begin()+j,make_pair(cluster(),minus[i]));
		  j=clusters.size();
		}
	    }
	  // 	  else
	  // 	    {
	  // 	      pushed=true;
	  // 	      clusters.push_back(make_pair(cluster(),minus[i]));
	  // 	      j=clusters.size();
	  // 	    }
	}
      if(!pushed)
	{
	  clusters.push_back(make_pair(cluster(),minus[i]));
	  //cerr << minus[i].positions.first << '\t' << minus[i].positions.second << '\n';
	}
    }
  /*
    for( unsigned i = 0 ; i < minus.size() ; ++i  )
    {
    for(unsigned j = 0 ; j < clusters.size() ; ++j )
    {
    if( clusters[j].first.positions.first != UMAX )
    {
    if( minus[i].positions.first < clusters[j].first.positions.first )
    {
    //add it
    clusters.insert(clusters.begin()+j,make_pair(cluster(),minus[i]));
    j=clusters.size();
    }
    else if ( clusters[j].second.positions.first != UMAX )
    {
    if( minus[i].positions.first <  clusters[j].second.positions.first )
    {
    clusters.insert(clusters.begin()+j,make_pair(cluster(),minus[i]));
    j=clusters.size();
    }
    }
    }
    else if( minus[i].positions.first < clusters[j].second.positions.first )
    {
    assert( clusters[j].second.positions.first != UMAX );
    clusters.insert(clusters.begin()+j,make_pair(cluster(),minus[i]));
    j=clusters.size();
		
    }
    else
    {
    clusters.push_back(make_pair(cluster(),minus[i]));
    j=clusters.size();
    }
    }
    }
  */
  /*
    for(unsigned i=0;i<clusters.size();++i)
    {
    cout << clusters[i].first.positions.first << '\t'
    << clusters[i].first.positions.second << '\t'
    << clusters[i].first.nreads << '\t'
    << clusters[i].second.positions.first << '\t'
    << clusters[i].second.positions.first << '\t'
    << clusters[i].second.nreads << '\n';
    }
  */
}

void output_results( //ofstream & out,
		    filtering_ostream & out,
		     const vector<pair<cluster,cluster> > & clusters, 
		     const string & chrom_label , 
		     const vector< pair<unsigned,unsigned> > & ref_te_chromo )
		     
{
  vector<pair<unsigned,unsigned> >::const_iterator mind;
  vector<pair<unsigned,unsigned> >::const_reverse_iterator mindr;
  unsigned mindist = numeric_limits<unsigned>::max();
  bool withinTE;
  out.flush();
  for(unsigned i=0;i<clusters.size();++i)
    {
      out.flush();
      out << chrom_label << '\t'
	  << clusters[i].first.nreads << '\t'
	  << clusters[i].second.nreads << '\t';
      if( clusters[i].first.positions.first == UMAX )
	{
	  out << "NA\t"
	      << "NA\t"
	      << "NA\t"
	      << "NA\t";
	}
      else
	{
	  out << clusters[i].first.positions.first << '\t'
	      << clusters[i].first.positions.second << '\t';
	  mind = find_if(ref_te_chromo.begin(),
			 ref_te_chromo.end(),
			 bind2nd(closest_plus(),clusters[i].first.positions.second));
	  if(mind != ref_te_chromo.end())
	    {
	      mindist = (mind->first < clusters[i].first.positions.first) ?
		clusters[i].first.positions.first-mind->first : 
		mind->first - clusters[i].first.positions.first;
	    }
	  withinTE = ( find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
			       bind2nd(within(),clusters[i].first.positions.first)) != ref_te_chromo.end() ||
		       find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
			       bind2nd(within(),clusters[i].first.positions.second)) 
		       != ref_te_chromo.end() );
	  if (mind != ref_te_chromo.end())
	    {
	      out << mindist << '\t';
	    }
	  else
	    {
	      out << "NA\t";
	    }
	  out << withinTE << '\t';
	}
      if( clusters[i].second.positions.first == UMAX )
	{
	  out << "NA\t"
	      << "NA\t"
	      << "NA\t"
	      << "NA" << endl;
	}
      else
	{
	  out << clusters[i].second.positions.first << '\t'
	      << clusters[i].second.positions.second << '\t';
	  mindr = find_if(ref_te_chromo.rbegin(),
			  ref_te_chromo.rend(),
			  bind2nd(closest_minus(),clusters[i].second.positions.first));
	  if(mindr != ref_te_chromo.rend())
	    {
	      mindist = (mindr->first < clusters[i].second.positions.second) ? 
		clusters[i].second.positions.second-mindr->first : mindr->first - clusters[i].second.positions.second;
	    }
	  withinTE = ( find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
			       bind2nd(within(),clusters[i].second.positions.second)) 
		       != ref_te_chromo.end() ||
		       find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
			       bind2nd(within(),clusters[i].second.positions.first))
		       != ref_te_chromo.end() );
	  if (mindr != ref_te_chromo.rend())
	    {
	      out << mindist << '\t';
	    }
	  else
	    {
	      out << "NA\t";
	    }
	  out << withinTE << endl;
	}
      out.flush();
    }
}

void reduce_ends( vector<cluster> & clusters,
		  const unsigned & INSERTSIZE )
{
  vector<cluster>::iterator i = clusters.end()-1,
    beg = clusters.begin(),j;

  //while(i>beg)
  while(i-clusters.begin() > 0)
    {
      bool merged = 0;
      assert(i>beg);
      for( j = i-1 ; !merged && j>=clusters.begin() ; --j )
	{
	  assert(j<i);
	  assert(j>=beg);
	  assert(beg==clusters.begin());
	  assert( i-beg > 0);
	  if( i->positions.first != UMAX && j->positions.first != UMAX )
	    {
	      assert(j-beg>=0);
	      if( (( max(i->positions.first,j->positions.first) -
		     min(i->positions.first,j->positions.first) ) <= INSERTSIZE ) ||
		  (( max(i->positions.first,j->positions.second) -
		     min(i->positions.first,j->positions.second) ) <= INSERTSIZE ) ||
		  (( max(i->positions.second,j->positions.second) -
		     min(i->positions.second,j->positions.second) ) <= INSERTSIZE ) ||
		  (( max(i->positions.second,j->positions.first) -
		     min(i->positions.second,j->positions.first) ) <= INSERTSIZE ) )
		{
		  j->positions.first = min(i->positions.first,
					   j->positions.first);
		  j->positions.second = max(i->positions.second,
					    j->positions.second);
		  j->nreads++;
		  clusters.erase(i);
		  i=clusters.end()-1;
		  j=i-1;
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
