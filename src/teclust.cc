/*
  Copyright 2010-2014 Kevin Thornton, University of California, Irvine
  
  This code is released under the terms of the GNU Public Licesne
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include <sstream>
#include <zlib.h>
#include <common.hpp>
#include <Sequence/IOhelp.hpp>

//Functions for this program
#include <teclust_objects.hpp>
#include <teclust_parseargs.hpp>
#include <teclust_scan_bamfile.hpp>
#include <teclust_phrapify.hpp>
#include <intermediateIO.hpp>

using namespace std;
using namespace Sequence;

using puu = pair<int32_t,int8_t>;

//const unsigned UMAX = std::numeric_limits<unsigned>::max();
const unsigned IMAX = std::numeric_limits<int32_t>::max();

refTEcont read_refdata( const teclust_params & p );
unordered_set<string> procUMM(const teclust_params & pars,
			      const refTEcont & reftes,
			      map<string,vector< puu > > * data);
void output_results(ostringstream & out,
		    const vector<pair<cluster,cluster> > & clusters, 
		    const string & chrom_label, 
		    const refTEcont & reftes);
void output_results_bedpe(ostringstream & out,
			  const vector<pair<cluster,cluster> > & clusters, 
			  const string & chrom_label, 
			  const refTEcont & reftes);
void cluster_data( vector<pair<cluster,cluster> > & clusters,
		   const vector<puu> & raw_data, 
		   const int32_t & INSERTSIZE, const int32_t & MDIST );
void reduce_ends( vector<cluster> & clusters,
		  const int32_t & INSERTSIZE );

int main( int argc, char ** argv )
{
  const teclust_params pars = teclust_parseargs(argc,argv);

  //Read in the locations of TEs in the reference
  auto refTEs = read_refdata(pars);
  /*
    Process the um_u and um_m files from the sample.  if refTEs is empty, parsedUMM contains the info for all U/M pairs.
    Otherwise, it contains only the info from U/M pairs where the M read hits a known TE in the reference.
  */
  //rawData = map {chromo x vector {start,strand}}
  map<string,vector< puu > > rawData;
  unordered_set<string> readPairs = procUMM(pars,refTEs,&rawData);
  /*
    Scan the BAM file to look for reads whose
    primary alignment hits a known TE in
    the reference, and whose mate is 
    mapped but does not hit a TE
  */
  scan_bamfile(pars,refTEs,&readPairs,&rawData);
  //Sort the raw data
  for( auto itr = rawData.begin();itr!=rawData.end();++itr )
    {
      sort(itr->second.begin(),itr->second.end(),
	   [](const puu & lhs, const puu & rhs) {
	     return lhs.first < rhs.first;
	   });
    }

  //Cluster the raw data and buffer results
  /*
  ostringstream out;
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
  */
  for( auto itr = rawData.begin() ; itr != rawData.end(); ++itr)
    {
      vector<pair<cluster,cluster> > clusters;
      cluster_data(clusters,itr->second,pars.INSERTSIZE,pars.MDIST);
      output_results_bedpe(out,clusters,
			   itr->first,
			   refTEs);
    }

  //write output
  gzFile gzout = gzopen(pars.outfile.c_str(),"w");
  if(gzout == NULL) 
    {
      cerr << "Error: "
	   << pars.outfile
	   << " could not be opened for writing.\n";
      exit(1);
    }
  if(! gzwrite(gzout,out.str().c_str(),unsigned(out.str().size())) )
    {
      cerr << "Error: gzwrite error at line " << __LINE__
	   << " of " << __FILE__ << '\n';
      exit(10);
    }
  gzclose(gzout);

  //Output the input to phrap, if desired
  phrapify( pars, out.str() );
}

refTEcont read_refdata( const teclust_params & p )
{
  refTEcont rv;
  
  if(p.reference_datafile.empty()) return rv;
  
  gzFile in = gzopen(p.reference_datafile.c_str(),"r");
  if(in == NULL )
    {
      cerr << "Error: could not open " 
	   << p.reference_datafile
	   << " for reading\n";
      exit(0);
    }

  do
    {
      auto line = IOhelp::gzreadline(in);
      if(!line.second) break;
      istringstream instream(line.first);
      string chrom;
      unsigned start,stop;
      instream >> chrom >> start >> stop >> ws;
      //The input file is BED, so start is 0 offset,
      //and stop is 1 offset.  Thus, we 
      //subtract 1 from stop to make both 0 offset
      rv[chrom].emplace_back( teinfo(start-1,stop) );
    }
  while(!gzeof(in));

  gzclose(in);

  //Sort the data
  for( auto __v = rv.begin() ; __v != rv.end() ; ++__v )
    {
      sort(__v->second.begin(),__v->second.end(),[](const teinfo & __l,const teinfo __r) {
	  return __l.start() < __r.start();
	});
    }
  return rv;
}


unordered_set<string> procUMM(const teclust_params & pars,
			      const refTEcont & reftes,
			      map<string,vector< puu > > * data)
{
  gzFile gzin = gzopen(pars.ummfile.c_str(),"r" );
  if(gzin == NULL)
    {
      cerr << "Error: "
	   << pars.ummfile
	   << " could not be opened for reading.\n";
      exit(1);
    }

  unordered_set<string> mTE; //"M" reads that map to a known TE in the refernce.  

  if (!reftes.empty() )
    {
      do
	{
	  auto name = gzreadCstr( gzin );
	  if(gzeof(gzin)) break;
	  if( name.second <= 0 ) 
	    {
	      cerr << "Error: gzread error at line "
		   << __LINE__ 
		   << " of " << __FILE__ << '\n';
	      exit(1);
	    }
	  auto chrom = gzreadCstr( gzin );
	  if( chrom.second <= 0 ) 
	    {
	      cerr << "Error: gzread error at line "
		   << __LINE__ 
		   << " of " << __FILE__ << '\n';
	      exit(1);
	    }
	  alnInfo alndata(gzin);
	  //Don't re-process a read if we already know it has a mapping to a TE
	  if( mTE.find(name.first) == mTE.end() )
	    {
	      auto __itr = reftes.find(chrom.first);
	      if( __itr != reftes.end() )
		{
		  //Default to the greedy algo of Cridland et al.
		  if( pars.greedy )
		    {
		      mTE.insert(name.first);
		    }
		  //Else, require that an M read overlap a TE
		  else if( find_if(__itr->second.cbegin(),
				   __itr->second.cend(),
				   [&](const teinfo & __t) {
				     return ( (alndata.start >= __t.start() && alndata.start <= __t.stop()) ||
					      (alndata.stop >= __t.start() && alndata.stop <= __t.stop()) );
				   }) != __itr->second.cend() )
		    {
		      //Then read hits a known TE
		      mTE.insert(name.first);
		    }
		}
	    }
	}
      while(!gzeof(gzin));
    }
  gzclose(gzin);

  //Now, get the Unique reads corresponding to TE-hitting M reads
  gzin = gzopen( pars.umufile.c_str(), "r" );
  if(gzin == NULL)
    {
      cerr << "Error: "
	   << pars.umufile
	   << " could not be opened for reading.\n";
      exit(1);
    }

  do
    {
      auto name = gzreadCstr( gzin );
      if(gzeof(gzin))break;
      if( name.second <= 0 ) 
	{
	  cerr << name.first << '\n';
	  cerr << "Error: gzread error at line "
	       << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      auto chrom = gzreadCstr( gzin );
      if( chrom.second <= 0 ) 
	{
	  cerr << "Error: gzread error at line "
	       << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      alnInfo alndata(gzin);
      if( reftes.empty() || (!reftes.empty() && mTE.find(name.first) != mTE.end()) )
	{
	  auto itr = data->find(chrom.first);
	  if( itr == data->end() )
	    {
	      data->insert(make_pair(chrom.first,vector<puu>(1,make_pair(alndata.start,alndata.strand))));
	    }
	  else
	    {
	      itr->second.push_back(make_pair(alndata.start,alndata.strand));
	    }
	}
    }
  while(!gzeof(gzin));
  gzclose(gzin);
  return mTE;
}


void cluster_data( vector<pair<cluster,cluster> > & clusters,
		   const vector<puu> & raw_data, 
		   const int32_t & INSERTSIZE, const int32_t & MDIST )
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
		  if( (max(raw_data[i].first,plus[j].positions.second)-
		       min(raw_data[i].first,plus[j].positions.second)) <= INSERTSIZE )
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
		  if( (max(raw_data[i].first,minus[j].positions.second)-
		       min(raw_data[i].first,minus[j].positions.second)) <= INSERTSIZE )
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
	    }
	}
    }
  
  reduce_ends( plus, INSERTSIZE );
  reduce_ends( minus, INSERTSIZE );

  auto close_enough_minus = [](const cluster & __minus,
			       const cluster & __plus,
			       const int32_t & __MDIST)
    {
      if( __minus.positions.first < __plus.positions.second ) return false;
      if( __minus.positions.first - __plus.positions.second <= __MDIST) return true;
      return false;
    };

  //now, we have to match up plus and minus based on MDIST
  vector<short> matched(plus.size(),0);
  for(unsigned i=0;i<plus.size();++i)
    {
      if(!matched[i])
	{
	  vector<cluster>::iterator j = find_if(minus.begin(),
						minus.end(),
						bind(close_enough_minus,placeholders::_1,plus[i],MDIST));
	  if( j != minus.end() )
	    {
	      //is there a better match in plus for this minus?
	      auto dist = j->positions.first-plus[i].positions.second;
	      unsigned winner = i;
	      for(unsigned k=i+1;k<plus.size();++k)
		{
		  if( close_enough_minus(*j,plus[k],MDIST) )
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

  //now, add in the minuses
  for( unsigned i = 0 ; i < minus.size() ; ++i  )
    {
      bool pushed = false;
      for(long unsigned j = 0 ; j < clusters.size() ; ++j )
	{
	  if( clusters[j].first.positions.first != IMAX )
	    {
	      if( minus[i].positions.first < clusters[j].first.positions.first )
		{
		  //add it
		  pushed=true;
		  clusters.insert(clusters.begin()+j,make_pair(cluster(),minus[i]));
		  j=clusters.size();
		}
	      else if ( clusters[j].second.positions.first != IMAX )
		{
		  if( minus[i].positions.first <  clusters[j].second.positions.first )
		    {
		      pushed=true;
		      clusters.insert(clusters.begin()+j,make_pair(cluster(),minus[i]));
		      j=clusters.size();
		    }
		}
	    }
	  else if( clusters[j].second.positions.first != IMAX )
	    {
	      assert( clusters[j].second.positions.first != IMAX );
	      if( minus[i].positions.first < clusters[j].second.positions.first )
		{
		  pushed=true;
		  clusters.insert(clusters.begin()+j,make_pair(cluster(),minus[i]));
		  j=clusters.size();
		}
	    }
	}
      if(!pushed)
	{
	  clusters.push_back(make_pair(cluster(),minus[i]));
	}
    }
}

void output_results( ostringstream & out,
		     const vector<pair<cluster,cluster> > & clusters, 
		     const string & chrom_label , 
		     const refTEcont & reftes )
{
  out.flush();

  auto closest_plus = [](const teinfo & __t,
			 const int32_t & rhs)
    {
      return  (__t.start() >= rhs || ( rhs >= __t.start() && rhs <= __t.stop() ) );
    };

  auto closest_minus = [](const teinfo & __t,
			  const int32_t & rhs)
    {
      return (__t.start() <= rhs || ( rhs >= __t.start() && rhs <= __t.stop() ) );
    };

  auto within = [](const teinfo & __t, const int32_t & start, const int32_t & stop)
    {
      bool A = start >= __t.start() && start <= __t.stop();
      bool B = stop>= __t.start() && stop <= __t.stop();
      return A||B;
    };

  auto refItr = reftes.find(chrom_label);
  for(unsigned i=0;i<clusters.size();++i)
    {
      out.flush();
      out << chrom_label << '\t'
	  << clusters[i].first.nreads << '\t'
	  << clusters[i].second.nreads << '\t';
      if( clusters[i].first.positions.first == IMAX )
	{
	  out << "-1\t"
	      << "-1\t"
	      << "-1\t"
	      << "-1\t";
	}
      else
	{
	  out << clusters[i].first.positions.first + 1 << '\t'
	      << clusters[i].first.positions.second + 1 << '\t';
	  int mindist = -1;
	  int withinTE = -1;
	  if(!reftes.empty())
	    {
	      auto mind = find_if( refItr->second.cbegin(), refItr->second.cend(),
				   bind(closest_plus,placeholders::_1,clusters[i].first.positions.second) );
	      if(mind != refItr->second.end())
		{
		  mindist = (mind->start() < clusters[i].first.positions.first) ?
		    clusters[i].first.positions.first-mind->start() : 
		    mind->start() - clusters[i].first.positions.first;
		}
	      auto __win = find_if(refItr->second.cbegin(),refItr->second.cend(),
				   bind(within,placeholders::_1,clusters[i].first.positions.first,clusters[i].first.positions.second));
	      withinTE = ( __win != refItr->second.cend() );
	    }
	  out << ((withinTE) ? 0 : mindist) << '\t' << withinTE << '\t';	  
	}
      if( clusters[i].second.positions.first == IMAX )
	{
	  out << "-1\t"
	      << "-1\t"
	      << "-1\t"
	      << "-1" << endl;
	}
      else
	{
	  out << clusters[i].second.positions.first + 1  << '\t'
	      << clusters[i].second.positions.second + 1 << '\t';
	  int mindist = -1;
	  int withinTE = -1;
	  if(!reftes.empty())
	    {
	      auto mindr = find_if(refItr->second.crbegin(),
				   refItr->second.crend(),
				   bind(closest_minus,placeholders::_1,clusters[i].second.positions.first));
	      if(mindr != refItr->second.crend())
		{
		  mindist = (mindr->start() < clusters[i].second.positions.second) ? 
		    clusters[i].second.positions.second-mindr->start() : mindr->start() - clusters[i].second.positions.second;
		}
	      auto __win = find_if(refItr->second.cbegin(),refItr->second.cend(),
				   bind(within,placeholders::_1,clusters[i].second.positions.first,clusters[i].second.positions.second));
	      withinTE = (__win != refItr->second.cend());
	    }
	  out << ((withinTE) ? 0 : mindist) << '\t' << withinTE << endl;
	}
      out.flush();
    }
}

void output_results_bedpe( ostringstream & out,
			   const vector<pair<cluster,cluster> > & clusters, 
			   const string & chrom_label , 
			   const refTEcont & reftes )
{
  out.flush();

  auto closest_plus = [](const teinfo & __t,
			 const int32_t & rhs)
    {
      return  (__t.start() >= rhs || ( rhs >= __t.start() && rhs <= __t.stop() ) );
    };

  auto closest_minus = [](const teinfo & __t,
			  const int32_t & rhs)
    {
      return (__t.start() <= rhs || ( rhs >= __t.start() && rhs <= __t.stop() ) );
    };

  auto within = [](const teinfo & __t, const int32_t & start, const int32_t & stop)
    {
      bool A = start >= __t.start() && start <= __t.stop();
      bool B = stop>= __t.start() && stop <= __t.stop();
      return A||B;
    };

  auto refItr = reftes.find(chrom_label);
  for(unsigned i=0;i<clusters.size();++i)
    {
      ostringstream xtra; //for the optional column
      // out.flush();
      // out << chrom_label << '\t'
      // 	  << clusters[i].first.nreads << '\t'
      // 	  << clusters[i].second.nreads << '\t';
      if( clusters[i].first.positions.first == IMAX )
	{
	  //Chrom, start, stop, unknown.
	  out << ".\t"
	      << "-1\t"
	      << "-1\t";
	}
      else
	{
	  out << chrom_label << '\t'
	      << clusters[i].first.positions.first << '\t'
	      << clusters[i].first.positions.second + 1 << '\t';
	  int mindist = -1;
	  int withinTE = -1;
	  if(!reftes.empty())
	    {
	      auto mind = find_if( refItr->second.cbegin(), refItr->second.cend(),
				   bind(closest_plus,placeholders::_1,clusters[i].first.positions.second) );
	      if(mind != refItr->second.end())
		{
		  mindist = (mind->start() < clusters[i].first.positions.first) ?
		    clusters[i].first.positions.first-mind->start() : 
		    mind->start() - clusters[i].first.positions.first;
		}
	      auto __win = find_if(refItr->second.cbegin(),refItr->second.cend(),
				   bind(within,placeholders::_1,clusters[i].first.positions.first,clusters[i].first.positions.second));
	      withinTE = ( __win != refItr->second.cend() );
	    }
	  xtra << ((withinTE) ? 0 : mindist) << '\t' << withinTE << '\t';
	  //out << ((withinTE) ? 0 : mindist) << '\t' << withinTE << '\t';	  
	}
      if( clusters[i].second.positions.first == IMAX )
	{
	  out << ".\t"
	      << "-1\t"
	      << "-1\t";
	}
      else
	{
	  out << chrom_label << '\t'
	      << clusters[i].second.positions.first << '\t'
	      << clusters[i].second.positions.second + 1 << '\t';
	  int mindist = -1;
	  int withinTE = -1;
	  if(!reftes.empty())
	    {
	      auto mindr = find_if(refItr->second.crbegin(),
				   refItr->second.crend(),
				   bind(closest_minus,placeholders::_1,clusters[i].second.positions.first));
	      if(mindr != refItr->second.crend())
		{
		  mindist = (mindr->start() < clusters[i].second.positions.second) ? 
		    clusters[i].second.positions.second-mindr->start() : mindr->start() - clusters[i].second.positions.second;
		}
	      auto __win = find_if(refItr->second.cbegin(),refItr->second.cend(),
				   bind(within,placeholders::_1,clusters[i].second.positions.first,clusters[i].second.positions.second));
	      withinTE = (__win != refItr->second.cend());
	    }
	  //out << ((withinTE) ? 0 : mindist) << '\t' << withinTE << endl;
	  xtra << ((withinTE) ? 0 : mindist) << '\t' << withinTE;// << endl;
	}
      out << "chr" << chrom_label << "_event_" << i << '\t'
	  <<  (clusters[i].first.nreads+clusters[i].second.nreads) << '\t'
	  << "+\t-\t"
	  << xtra.str() << '\n';
      out.flush();
    }
}

void reduce_ends( vector<cluster> & clusters,
		  const int32_t & INSERTSIZE )
{
  vector<cluster>::iterator i = clusters.end()-1,
    beg = clusters.begin(),j;

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
	  if( i->positions.first != IMAX && j->positions.first != IMAX )
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


