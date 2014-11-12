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
#include <unordered_set>
#include <vector>
#include <functional>
#include <algorithm>
#include <limits>
#include <cassert>
#include <sstream>
#include <zlib.h>
#include <Sequence/IOhelp.hpp>
#include <Sequence/bamreader.hpp>
#include <Sequence/bamrecord.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace Sequence;
using namespace boost::program_options;

// DEFINITIONS OF DATA TYPES
using puu = pair<unsigned,unsigned>;

const unsigned UMAX = std::numeric_limits<unsigned>::max();

// struct teinfo
// //! Where TEs are in the reference genome
// {
//   string chrom;
//   unsigned start,stop;
//   teinfo( const string & c, const unsigned & st, const unsigned & sto) :
//   chrom(c),start(st),stop(sto)
//   {
//   }
// };

struct teinfo : public pair<unsigned,unsigned>
{
  unsigned start() const { return this->first; }
  unsigned stop() const { return this->second; }
  teinfo( unsigned __s, unsigned __st ) : pair<unsigned,unsigned>(move(__s),move(__st))
  {
  }
};

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

struct params //Command-line parameter options
{
  /*
    The reference TE positions, 
    output file name, 
    input bam file name, 
    file name listing the FASTQ files
    the um_u file that is the output of cluster_cnv for the sample
    the um_m files that is the output of cluster_cnv for the sample
    The output of this program run on the reference genome, if available
    */
  string  reference_datafile, outfile, bamfile, readfile, umufile, ummfile;//,teclust_ref;
  /*
    Upper limit on insert size distribution
    Maximum distance used for matching up left and right ends of putative TE calls
   */
  unsigned INSERTSIZE,MDIST;
  params();
}; 

params::params() : reference_datafile(string()),
		   outfile(string()),
		   bamfile(string()),
		   readfile(string()),
		   umufile(string()),
		   ummfile(string()),
		   INSERTSIZE(UMAX),
		   MDIST(UMAX)
{
}

using refTEcont = map<string,vector<teinfo> >;
//DEFINITION OF FUNCTIONS
params parseargs(const int argc, char ** argv);
vector<pair<string,int32_t> > make_lookup(const bamreader & reader);
refTEcont read_refdata( const params & p );
void procUMM(const refTEcont & reftes,
	     const string & umufilename,
	     const string & ummfilename,
	     map<string,vector< puu > > * data);
void output_results(ostringstream & out,
		    const vector<pair<cluster,cluster> > & clusters, 
		    const string & chrom_label, 
		    const refTEcont & reftes);
		    //const vector< teinfo > & reftes );
void scan_bamfile(const params & p,
		  const refTEcont & refTEs,
		  map<string,vector< puu > > * data);
//OLD

void cluster_data( vector<pair<cluster,cluster> > & clusters,
		   const vector<puu> & raw_data, 
		   const unsigned & INSERTSIZE, const unsigned & MDIST );
void reduce_ends( vector<cluster> & clusters,
		  const unsigned & INSERTSIZE );
// void output_results(ostringstream & out,
// 		    const vector<pair<cluster,cluster> > & clusters, 
// 		    const string & chrom_label, 
// 		    const vector< pair<unsigned,unsigned> > & ref_te_chromo);

// void read_raw_data(gzFile gzin,
// 		   map<unsigned,vector<puu> > & raw_data,
// 		   vector<pair<string,unsigned> > * chrom_labels )
// /*
//   input data are read from STDIN, and look like:
//   position chrom strand
  
// */
// {
//   unsigned pos,chrom,strand,dummy=0;
//   string chrom_label;
//   map<unsigned,vector<puu> >::iterator itr;
//   unsigned nread=0;

//   do
//     {
//       auto nextline = Sequence::IOhelp::gzreadline(gzin);
//       istringstream in(nextline.first);
//       in >> pos >> chrom_label >> strand  >> ws;
//       ++nread;
//       chrom = update_lookup(chrom_labels,&dummy,chrom_label);
//       itr = raw_data.find(chrom);

//       if( itr == raw_data.end() ) //new chromosome
// 	{
// 	  raw_data.insert( make_pair(chrom,
// 				     vector<puu>(1,puu(pos,strand))) );
// 	}
//       else //add to existing data for chromosome
// 	{
// 	  itr->second.push_back(puu(pos,strand));
// 	}
//     } while(!gzeof(gzin));
//   cerr << nread << " lines processed\n";
//   //sort data per chromosome by position
// }

int main( int argc, char ** argv )
{
  const params pars = parseargs(argc,argv);

  //Read in the locations of TEs in the reference
  auto refTEs = read_refdata(pars);

  /*
    Process the um_u and um_m files from the sample.  if refTEs is empty, parsedUMM contains the info for all U/M pairs.
    Otherwise, it contains only the info from U/M pairs where the M read hits a known TE in the reference.
  */
  map<string,vector< pair<unsigned,unsigned> > > rawData;
  procUMM(refTEs,pars.umufile,pars.ummfile,&rawData);
  
  scan_bamfile(pars,refTEs,&rawData);


  //Sort the raw data
  for( auto itr = rawData.begin();itr!=rawData.end();++itr )
    {
      sort(itr->second.begin(),itr->second.end(),
	   [](const puu & lhs, const puu & rhs) {
	     return lhs.first < rhs.first;
	   });
    }

  //Cluster the raw data and buffer results
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
  for( auto itr = rawData.begin() ; itr != rawData.end(); ++itr)
    {
      vector<pair<cluster,cluster> > clusters;
      cluster_data(clusters,itr->second,pars.INSERTSIZE,pars.MDIST);
      output_results(out,clusters,
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
  if(! gzwrite(gzout,out.str().c_str(),out.str().size()) )
    {
      cerr << "Error: gzwrite error at line " << __LINE__
	   << " of " << __FILE__ << '\n';
      exit(10);
    }
  gzclose(gzout);
}

params parseargs(const int argc, char ** argv)
{
  params rv;
  options_description desc("Cluster reads into putative transposable element calls.\nUsage: tclust -h to see help");
  desc.add_options()
    ("help,h", "Produce help message")
    ("bamfile,b",value<string>(&rv.bamfile),"BAM file name (optional)")
    ("tepos,t",value<string>(&rv.reference_datafile),"File containing positions of TEs in reference genome (optional)")
    ("outfile,o",value<string>(&rv.outfile),"Output file name for clusters (required)")
    ("umu,u",value<string>(&rv.umufile),"The um_u output file for the sample generated by process_readmappings (required)")
    ("umm,m",value<string>(&rv.ummfile),"The um_u output file for the sample generated by process_readmappings (required)")
    //("umuref,r",value<string>(&rv.teclust_ref),"The output of this program run on a reference genome, if available (optional)")
    ("isize,i",value<unsigned>(&rv.INSERTSIZE),"Upper limit on insert size distribution, e.g. from bwa_mapdistance. (required)")
    ("mdist,M",value<unsigned>(&rv.MDIST)->default_value(1000u),"Max. distance for joining up left and right ends of putative TEs (required)")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if( argc == 1 || 
      vm.count("help") ||
      !vm.count("outfile") ||
      !vm.count("umu") ||
      !vm.count("umm") ||
      !vm.count("isize") ||
      !vm.count("mdist") )
    {
      cerr << desc << '\n';
      exit(0);
    }

  return rv;
}

vector<pair<string,int32_t> > make_lookup(const bamreader & reader)
{
  vector<pair<string,int32_t> > rv;
  std::copy(reader.ref_cbegin(),reader.ref_cend(),
	    back_inserter(rv));
  return rv;
}

refTEcont read_refdata( const params & p )
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
      //We subtract 1 to convert to a 0-offest system
      rv[chrom].emplace_back( teinfo(start-1,stop-1) );
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
  // for_each(rv.begin(),rv.end(),
  // 	   []( const pair<string,vector<teinfo>> & data )
  // 	   {
  // 	     for_each(data.second.begin(),data.second.end(),[&](const teinfo & __t)
  // 		      {
  // 			cout << data.first << ' ' << __t.start() << ' ' << __t.stop() << '\n';
  // 		      });
  // 	   });
  // exit(0);
  return rv;
}


void procUMM(const refTEcont & reftes,
	     const string & umufilename,
	     const string & ummfilename,
	     map<string,vector<pair<unsigned,unsigned> > > * data)
{
  gzFile gzin = gzopen(ummfilename.c_str(),"r" );
  if(gzin == NULL)
    {
      cerr << "Error: "
	   << ummfilename
	   << " could not be opened for reading.\n";
      exit(1);
    }

  unordered_set<string> mTE; //"M" reads that map to a known TE in the refernce.  

  if (!reftes.empty() )
    {
      do
	{
	  auto line = IOhelp::gzreadline(gzin);
	  if(!line.second) break;
	  istringstream in(line.first);
	  string name,chrom;
	  unsigned mapq,start,stop;
	  //Only parse what we need to
	  in >> name >> mapq >> chrom >> start >> stop >> ws;
	  
	  //Don't re-process a read if we already know it has a mapping to a TE
	  if( mTE.find(name) != mTE.end() )
	    {
	      auto __itr = reftes.find(chrom);
	      if( find_if(__itr->second.cbegin(),
			  __itr->second.cend(),
			  [&](const teinfo & __t) {
			    return ( (start >= __t.start() && start <= __t.stop()) ||
				     (stop >= __t.start() && stop <= __t.stop()) );
			  }) != __itr->second.cend() )
		{
		  //Then read hits a known TE
		  mTE.insert(name);
		}
	    }
	}
      while(!gzeof(gzin));
    }

  gzclose(gzin);

  //Now, get the Unique reads corresponding to TE-hitting M reads
  gzin = gzopen( umufilename.c_str(), "r" );
  if(gzin == NULL)
    {
      cerr << "Error: "
	   << umufilename
	   << " could not be opened for reading.\n";
      exit(1);
    }

  do
    {
      auto line = IOhelp::gzreadline(gzin);
      if(!line.second) break;
      istringstream in(line.first);
      string name,chrom;
      unsigned mapq,start,stop,strand;
      //Only parse what we need to
      in >> name >> mapq >> chrom >> start >> stop >> strand >> ws;
      if( reftes.empty() || (!reftes.empty() && mTE.find(name) != mTE.end()) )
	{
	  auto itr = data->find(chrom);
	  if( itr == data->end() )
	    {
	      data->insert(make_pair(chrom,vector<pair<unsigned,unsigned> >(1,make_pair(start,strand))));
	    }
	  else
	    {
	      itr->second.push_back(make_pair(start,strand));
	    }
	}
    }
  while(!gzeof(gzin));

  gzclose(gzin);
}

void scan_bamfile(const params & p,
		  const refTEcont & refTEs,
		  map<string,vector< puu > > * data)
{
  if( refTEs.empty() || p.bamfile.empty() ) return;

  bamreader reader(p.bamfile.c_str());
  if(! reader )
    {
      cerr << "Error: " 
	   << p.bamfile
	   << " could not be opened for reading.\n";
      exit(0);
    }

  auto lookup = make_lookup(reader);
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
			       const unsigned & __MDIST)
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
	      unsigned dist = j->positions.first-plus[i].positions.second;
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
		     //const vector< teinfo > & reftes )
		     
{
  vector<pair<unsigned,unsigned> >::const_iterator mind;
  vector<pair<unsigned,unsigned> >::const_reverse_iterator mindr;
  int32_t mindist = -1;
  int withinTE = -1;
  out.flush();

  auto closest_plus = [](const teinfo & __t,
			 const unsigned & rhs)
    {
      return  (__t.start() >= rhs || ( rhs >= __t.start() && rhs <= __t.stop() ) );
    };

  auto closest_minus = [](const teinfo & __t,
			  const unsigned & rhs)
    {
      return (__t.start() <= rhs || ( rhs >= __t.start() && rhs <= __t.stop() ) );
    };

  auto refItr = reftes.find(chrom_label);
  for(unsigned i=0;i<clusters.size();++i)
    {
      out.flush();
      out << chrom_label << '\t'
	  << clusters[i].first.nreads << '\t'
	  << clusters[i].second.nreads << '\t';
      if( clusters[i].first.positions.first == UMAX )
	{
	  out << "-1\t"
	      << "-1\t"
	      << "-1\t"
	      << "-1\t";
	}
      else
	{
	  out << clusters[i].first.positions.first << '\t'
	      << clusters[i].first.positions.second << '\t';
	  // mind = find_if(ref_te_chromo.begin(),
	  // 		 ref_te_chromo.end(),
	  // 		 [&]( const pair<unsigned,unsigned> & lhs ) {
	  // 		   //This is the old closest_plus function object from 0.1.0
	  // 		   //Finds the closest TE in the reference 3' of this position
	  // 		   const unsigned rhs = clusters[i].first.positions.second;
	  // 		   return ( lhs.first >= rhs || ( rhs >= lhs.first && rhs <= lhs.second ) );
	  // 		 });
	  mindist = -1;
	  withinTE = -1;
	  if(!reftes.empty())
	    {
	      auto mind = find_if( refItr->second.cbegin(), refItr->second.cend(),
				   bind(closest_plus,placeholders::_1,clusters[i].first.positions.second) );
	      // [&](const teinfo & __t) {
	      // 	 const unsigned rhs = clusters[i].first.positions.second;
	      // 	 return (chrom_label == __t.chrom && ( __t.start >= rhs || ( rhs >= __t.start && rhs <= __t.stop ) ));
	      // });
	      if(mind != refItr->second.end())//ref_te_chromo.end())
		{
		  mindist = (mind->start() < clusters[i].first.positions.first) ?
		    clusters[i].first.positions.first-mind->start() : 
		    mind->start() - clusters[i].first.positions.first;
		}
	    
	  // withinTE = ( find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
	  // 		       [&](const pair<unsigned,unsigned> & refTE) {
	  // 			 return clusters[i].first.positions.first >= refTE.first ||
	  // 			 clusters[i].first.positions.first <= refTE.second;
	  // 		       }) != ref_te_chromo.end() ||
	  // 		       //bind2nd(within(),clusters[i].first.positions.first)) != ref_te_chromo.end() ||
	  // 	       find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
	  // 		       [&](const pair<unsigned,unsigned> & refTE) {
	  // 			 return clusters[i].first.positions.second >= refTE.first || 
	  // 			 clusters[i].first.positions.second <= refTE.second;
	  // 		       })
	  // 		       //bind2nd(within(),clusters[i].first.positions.second)) 
	  // 	       != ref_te_chromo.end() );
	  withinTE = reftes.empty() ? false : ( find_if(refItr->second.cbegin(),refItr->second.cend(),
							[&](const teinfo & __t) {
							  bool B = (clusters[i].first.positions.first >= __t.start() || clusters[i].first.positions.first <= __t.stop());
							  bool C = (clusters[i].first.positions.second >= __t.start() || clusters[i].first.positions.second <= __t.stop());
							  return (B||C);
							  /*
							    return __t.chrom == chrom_label &&
							    ((clusters[i].first.positions.first >= __t.start ||
							    clusters[i].first.positions.first <= __t.stop) ||
							    (clusters[i].first.positions.second >= __t.start ||
							    clusters[i].first.positions.second <= __t.stop));
							  */
							}) != refItr->second.cend() );
	    }
	  out << mindist << '\t' << withinTE << '\t';	  
	  // if (mind != refItr->second.end())//ref_te_chromo.end())
	  //   {
	  //     out << mindist << '\t';
	  //   }
	  // else
	  //   {
	  //     out << "-1\t";
	  //   }
	  //out << withinTE << '\t';
	  //out << ((!reftes.empty())?int(withinTE):-1) << '\t';
	}
      if( clusters[i].second.positions.first == UMAX )
	{
	  out << "-1\t"
	      << "-1\t"
	      << "-1\t"
	      << "-1" << endl;
	}
      else
	{
	  out << clusters[i].second.positions.first << '\t'
	      << clusters[i].second.positions.second << '\t';
	  /*
	  mindr = find_if(ref_te_chromo.rbegin(),
			  ref_te_chromo.rend(),
			  //This is the old closest_minus from 0.1.0
			  //Finds the closest reference TE 5' of this position
			  [&](const pair<unsigned,unsigned> & lhs) {
			    const unsigned rhs = clusters[i].second.positions.first;
			    return ( lhs.second <= rhs || ( rhs >= lhs.first && rhs <= lhs.second ) );
			  });
	  */
	  mindist = -1;
	  withinTE = -1;
	  if(!reftes.empty())
	    {
	      auto mindr = find_if(refItr->second.crbegin(),
				   refItr->second.crend(),
				   bind(closest_minus,placeholders::_1,clusters[i].second.positions.first));
	      /*
			       [&](const teinfo __t){
				 const unsigned rhs = clusters[i].second.positions.first;
				 return chrom_label == __t.chrom &&
				 (__t.stop <= rhs || ( rhs >= __t.start && rhs <= __t.stop ));
			       });
	      */
	      if(mindr != refItr->second.crend())//ref_te_chromo.rend())
		{
		  mindist = (mindr->start() < clusters[i].second.positions.second) ? 
		    clusters[i].second.positions.second-mindr->start() : mindr->start() - clusters[i].second.positions.second;
		}
	  // withinTE = ( find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
	  // 		       [&](const pair<unsigned,unsigned> & refTE) {
	  // 			 return clusters[i].second.positions.second >= refTE.first ||
	  // 			 clusters[i].second.positions.second <= refTE.second;
	  // 		       })
	  // 		       //bind2nd(within(),clusters[i].second.positions.second)) 
	  // 	       != ref_te_chromo.end() ||
	  // 	       find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
	  // 		       [&](const pair<unsigned,unsigned> & refTE) {
	  // 			 return clusters[i].second.positions.first >= refTE.first ||
	  // 			 clusters[i].second.positions.first <= refTE.second;
	  // 		       })
	  // 		       //bind2nd(within(),clusters[i].second.positions.first))
	  // 	       != ref_te_chromo.end() );
	      withinTE = (find_if(refItr->second.cbegin(),refItr->second.cend(),
				  [&](const teinfo & __t) {
				    //bool A = chrom_label == __t.chrom;
				    bool B = (clusters[i].first.positions.first >= __t.start() || clusters[i].first.positions.first <= __t.stop());
				    bool C = (clusters[i].first.positions.second >= __t.start() || clusters[i].first.positions.second <= __t.stop());
				    return (B||C);
				    // return chrom_label == __t.chrom &&
				    // ((clusters[i].first.positions.first >= __t.start ||
				    //   clusters[i].first.positions.first <= __t.stop) ||
				    //  (clusters[i].first.positions.second >= __t.start ||
				    //   clusters[i].first.positions.second <= __t.stop));
				  }) != refItr->second.cend());
	    }
	  out << mindist << '\t' << withinTE << endl;
	  // if (mindr != refItr->second.crend()) // ref_te_chromo.rend())
	  //   {
	  //     out << mindist << '\t';
	  //   }
	  // else
	  //   {
	  //     out << "-1\t";
	  //   }
	  // out << ((!reftes.empty())?int(withinTE):-1) << endl;
	}
      out.flush();
    }
}

//OLD
// void output_results( ostringstream & out,
// 		     const vector<pair<cluster,cluster> > & clusters, 
// 		     const string & chrom_label , 
// 		     const vector< pair<unsigned,unsigned> > & ref_te_chromo )
		     
// {
//   vector<pair<unsigned,unsigned> >::const_iterator mind;
//   vector<pair<unsigned,unsigned> >::const_reverse_iterator mindr;
//   unsigned mindist = numeric_limits<unsigned>::max();
//   bool withinTE;
//   out.flush();
//   for(unsigned i=0;i<clusters.size();++i)
//     {
//       out.flush();
//       out << chrom_label << '\t'
// 	  << clusters[i].first.nreads << '\t'
// 	  << clusters[i].second.nreads << '\t';
//       if( clusters[i].first.positions.first == UMAX )
// 	{
// 	  out << "NA\t"
// 	      << "NA\t"
// 	      << "NA\t"
// 	      << "NA\t";
// 	}
//       else
// 	{
// 	  out << clusters[i].first.positions.first << '\t'
// 	      << clusters[i].first.positions.second << '\t';
// 	  mind = find_if(ref_te_chromo.begin(),
// 			 ref_te_chromo.end(),
// 			 [&]( const pair<unsigned,unsigned> & lhs ) {
// 			   //This is the old closest_plus function object from 0.1.0
// 			   //Finds the closest TE in the reference 3' of this position
// 			   const unsigned rhs = clusters[i].first.positions.second;
// 			   return ( lhs.first >= rhs || ( rhs >= lhs.first && rhs <= lhs.second ) );
// 			 });
// 	  if(mind != ref_te_chromo.end())
// 	    {
// 	      mindist = (mind->first < clusters[i].first.positions.first) ?
// 		clusters[i].first.positions.first-mind->first : 
// 		mind->first - clusters[i].first.positions.first;
// 	    }
// 	  withinTE = ( find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
// 			       [&](const pair<unsigned,unsigned> & refTE) {
// 				 return clusters[i].first.positions.first >= refTE.first ||
// 				 clusters[i].first.positions.first <= refTE.second;
// 			       }) != ref_te_chromo.end() ||
// 			       //bind2nd(within(),clusters[i].first.positions.first)) != ref_te_chromo.end() ||
// 		       find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
// 			       [&](const pair<unsigned,unsigned> & refTE) {
// 				 return clusters[i].first.positions.second >= refTE.first || 
// 				 clusters[i].first.positions.second <= refTE.second;
// 			       })
// 			       //bind2nd(within(),clusters[i].first.positions.second)) 
// 		       != ref_te_chromo.end() );
// 	  if (mind != ref_te_chromo.end())
// 	    {
// 	      out << mindist << '\t';
// 	    }
// 	  else
// 	    {
// 	      out << "NA\t";
// 	    }
// 	  out << withinTE << '\t';
// 	}
//       if( clusters[i].second.positions.first == UMAX )
// 	{
// 	  out << "NA\t"
// 	      << "NA\t"
// 	      << "NA\t"
// 	      << "NA" << endl;
// 	}
//       else
// 	{
// 	  out << clusters[i].second.positions.first << '\t'
// 	      << clusters[i].second.positions.second << '\t';
// 	  mindr = find_if(ref_te_chromo.rbegin(),
// 			  ref_te_chromo.rend(),
// 			  //This is the old closest_minus from 0.1.0
// 			  //Finds the closest reference TE 5' of this position
// 			  [&](const pair<unsigned,unsigned> & lhs) {
// 			    const unsigned rhs = clusters[i].second.positions.first;
// 			    return ( lhs.second <= rhs || ( rhs >= lhs.first && rhs <= lhs.second ) );
// 			  });
// 	  if(mindr != ref_te_chromo.rend())
// 	    {
// 	      mindist = (mindr->first < clusters[i].second.positions.second) ? 
// 		clusters[i].second.positions.second-mindr->first : mindr->first - clusters[i].second.positions.second;
// 	    }
// 	  withinTE = ( find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
// 			       [&](const pair<unsigned,unsigned> & refTE) {
// 				 return clusters[i].second.positions.second >= refTE.first ||
// 				 clusters[i].second.positions.second <= refTE.second;
// 			       })
// 			       //bind2nd(within(),clusters[i].second.positions.second)) 
// 		       != ref_te_chromo.end() ||
// 		       find_if(ref_te_chromo.begin(),ref_te_chromo.end(),
// 			       [&](const pair<unsigned,unsigned> & refTE) {
// 				 return clusters[i].second.positions.first >= refTE.first ||
// 				 clusters[i].second.positions.first <= refTE.second;
// 			       })
// 			       //bind2nd(within(),clusters[i].second.positions.first))
// 		       != ref_te_chromo.end() );
// 	  if (mindr != ref_te_chromo.rend())
// 	    {
// 	      out << mindist << '\t';
// 	    }
// 	  else
// 	    {
// 	      out << "NA\t";
// 	    }
// 	  out << withinTE << endl;
// 	}
//       out.flush();
//     }
// }

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
