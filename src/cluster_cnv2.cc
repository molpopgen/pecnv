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

#include <string>
#include <sstream>
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <zlib.h>
#include <intermediateIO.hpp>
#include <boost/program_options.hpp>
#include <sys/stat.h>

using namespace std;
using namespace boost::program_options;

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
	     const short & _strand1,
	     const short & _strand2) : a(__a),
				       aS(__aS),
				       b(__b),
				       bS(__bS),
				       readname( __readname ),
				       strand1(_strand1),strand2(_strand2)
  {
  }
};

auto order_clusters = []( const vector<vector<linkeddata>::const_iterator> & a,
			  const vector<vector<linkeddata>::const_iterator> & b ) 
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

using cluster_container = vector< vector<vector<linkeddata>::const_iterator> >;
using lvector = vector<linkeddata>;
using putCNVs = map<string,lvector>;

cluster_container cluster_linked( const lvector & raw,
				  const unsigned & mdist );

void reduce_clusters( cluster_container & clusters,
		      const unsigned & mdist );

unsigned mindist(const unsigned & st1,
		 const unsigned & stp1,
		 const unsigned & st2,
		 const unsigned & stp2);

bool pair_should_cluster( lvector::const_iterator & pair,
			  vector<lvector::const_iterator> & cluster,
			  const unsigned & mdist);

bool unique_positions(const lvector & data,
		      const unsigned & start,
		      const unsigned & start2);

void write_clusters( gzFile o,
		     const string & chrom1,
		     const string & chrom2,
		     const cluster_container & clusters,
		     unsigned * eventid );

void write_clusters_bedpe( gzFile o,
			   const string & sampleID,
			   const string & eventType,
			   const string & chrom1,
			   const string & chrom2,
			   const cluster_container & clusters,
			   unsigned * eventid );

void read_data(putCNVs & raw_div,
	       putCNVs & raw_par,
	       map<string,putCNVs > & raw_ul,
	       const char * filename,
	       const int8_t & min_mqual,
	       const int16_t & max_mm,
	       const int16_t & max_gap);

struct cluster_cnv_params
{
  string sampleID;
  int8_t min_mqual;
  int16_t max_mm,max_gap;
  unsigned mdist;
  string divfile,parfile,ulfile;
  vector<string> infiles;
};

cluster_cnv_params clusterCNV_parseargs(int argc, char ** argv);

int cluster_cnv_main(int argc, char ** argv)
{
  auto pars = clusterCNV_parseargs(argc, argv);
  // int argn=1;
  // if( argc < 9 )
  //   {
  //     cerr << "usage: "
  // 	   << argv[0]
  // 	   << " sampleID min_quality max_mm max_gap insert_size outfile_div outfile_par outfile_ul "
  // 	   << "structural_file1a structural_file1b ... structural_fileNa structural_fileNb\n";
  //     exit(0);
  //   }
  // const string sampleID(argv[argn++]);
  // const int8_t min_mqual = stoi(argv[argn++]);
  // const int16_t max_mm = stoi(argv[argn++]);
  // const int16_t max_gap = stoi(argv[argn++]);
  // const unsigned mdist = stoi(argv[argn++]);
  // const char * divfile = argv[argn++];
  // const char * parfile = argv[argn++];
  // const char * ulfile = argv[argn++];

  //make sure output files are writable
  //const string header = "id\tchrom1\tcoverage\tstrand1\tstart1\tstop1\tchrom2\tstrand2\tstart2\tstop2\treads";

  gzFile divstream = gzopen(pars.divfile.c_str(),"wb");
  if(divstream==NULL) {
    cerr << "Error: could not open "
	 << pars.divfile
	 << " for writing\n";
  }
  /*
 if( gzprintf(divstream,"%s\n",header.c_str()) <= 0 )
    {
      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }
  */
  gzFile parstream = gzopen(pars.parfile.c_str(),"wb");
  if(parstream == NULL) {
    cerr << "Error: could not open "
	 << pars.parfile << " for writing\n";
    exit(1);
  }
  /*
  if (gzprintf(parstream,"%s\n",header.c_str()) <= 0 )
    {
      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }
  */
  gzFile ulstream = gzopen(pars.ulfile.c_str(),"wb");
  if(ulstream == NULL)
    {
      cerr << "Error: could not open "
	   << pars.ulfile
	   << " for writing\n";
      exit(1);
    }
  /*
  if (gzprintf(ulstream,"%s\n",header.c_str()) <= 0 )
    {
      cerr << "Error: gzprintf error encountered at line " << __LINE__ 
	   << " of " << __FILE__ << '\n';
      exit(1);
    }
  */
  map<string, lvector > raw_div;
  map<string, lvector > raw_par;
  map<string, putCNVs > raw_ul;
  //for(int i = argn;i<argc;++i)//i+=2)
  for(unsigned i = 0 ; i < pars.infiles.size() ; ++i )
    {
      cerr << "processing " << argv[i] << '\n';
      read_data(raw_div,raw_par,raw_ul,
		pars.infiles[i].c_str(),
		pars.min_mqual,
		pars.max_mm,
		pars.max_gap);
    }

  unsigned eventid=0;
  cerr << "clustering div\n";
  for(putCNVs::iterator itr = raw_div.begin();
      itr != raw_div.end();++itr)
    {
      sort(itr->second.begin(),
	   itr->second.end(),
	   [](const linkeddata & lhs, const linkeddata & rhs){
	     return lhs.a < rhs.a && lhs.b < rhs.b;
	   });
      cluster_container clusters = cluster_linked(itr->second,pars.mdist);
      cluster_container clusters2(clusters);
      sort(clusters.begin(),clusters.end(),order_clusters);
      write_clusters_bedpe( divstream, 
			    pars.sampleID,
			    string("div"),
			    itr->first,itr->first,
			    clusters,&eventid );
    }

  cerr << "clustering par\n";
  eventid=0;
  for(putCNVs::iterator itr = raw_par.begin();
      itr != raw_par.end();++itr)
    {
      sort(itr->second.begin(),
	   itr->second.end(),
	   [](const linkeddata & lhs, const linkeddata & rhs){
	     return lhs.a < rhs.a && lhs.b < rhs.b;
	   });
      cluster_container clusters = cluster_linked(itr->second,pars.mdist);
      sort(clusters.begin(),clusters.end(),order_clusters);
      write_clusters_bedpe( parstream, 
			    pars.sampleID,
			    string("par"),
			    itr->first,
			    itr->first,
			    clusters,&eventid );
    }

  cerr << "clustering ul\n";
  eventid=0;
  for( map<string, putCNVs >::iterator itr = raw_ul.begin() ;
       itr != raw_ul.end() ; ++itr )
    {
      for( putCNVs::iterator itr2 = itr->second.begin() ; 
	   itr2 != itr->second.end() ; ++itr2 )
	{
	  assert(itr->first < itr2->first);
	  sort(itr2->second.begin(),itr2->second.end(),
	       [](const linkeddata & lhs, const linkeddata & rhs){
		 return lhs.a < rhs.a && lhs.b < rhs.b;
	       });
	  cluster_container clusters = cluster_linked(itr2->second,pars.mdist);
	  sort(clusters.begin(),clusters.end(),order_clusters);
	  write_clusters_bedpe( ulstream, 
				pars.sampleID,
				string("unl"),
				itr->first,
				itr2->first,
				clusters,&eventid );
	}
    }
  gzclose(parstream);
  gzclose(ulstream);
  gzclose(divstream);
}

cluster_cnv_params clusterCNV_parseargs(int argc, char ** argv)
{
  cluster_cnv_params rv;

  // const string sampleID(argv[argn++]);
  // const int8_t min_mqual = stoi(argv[argn++]);
  // const int16_t max_mm = stoi(argv[argn++]);
  // const int16_t max_gap = stoi(argv[argn++]);
  // const unsigned mdist = stoi(argv[argn++]);
  // const char * divfile = argv[argn++];
  // const char * parfile = argv[argn++];
  // const char * ulfile = argv[argn++];
  int mqual;
  options_description desc("pecnv cnvclust: cluster divergent, parallel, and unlinked read pairs into putative CNV calls");
  desc.add_options()
    ("help,h", "Produce help message")
    ("infiles,i",value<vector<string> >()->multitoken(),"Input files.  The input files are the output from the pecnv process subcommand")
    ("mqual,m",value<int>(&mqual)->default_value(30),"Minimum mapping quality for a read to be included")
    ("maxdist,d",value<unsigned>(&rv.mdist),"Upper limit of insert size distribution")
    ("sample,s",value<string>(&rv.sampleID)->default_value("sample"),"Unique label/name for the sample")
    ("maxmm,M",value<int16_t>(&rv.max_mm)->default_value(3),"Max no. mismatches in a read for it to be included")
    ("maxgap,g",value<int16_t>(&rv.max_gap)->default_value(0),"Max no. alignment gaps in a read for it to be included")
    ("divfile,D",value<string>(&rv.divfile)->default_value("div_clusters.gz"),"Output file for divergent clusters")
    ("parfile,P",value<string>(&rv.parfile)->default_value("par_clusters.gz"),"Output file for parallel clusters")
    ("parfile,U",value<string>(&rv.ulfile)->default_value("unl_clusters.gz"),"Output file for unlinked clusters")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);


  if( argc == 1 || 
      vm.count("help") ||
      !vm.count("infiles") ||
      !vm.count("mqual") ||
      !vm.count("maxdist") )
    {
      cerr << desc << '\n';
      exit(1);
    }
  if(!vm.count("mqual"))
    {
      rv.min_mqual = int8_t(mqual);
    }

  rv.infiles = vm["infiles"].as<vector<string> >();

  for(unsigned i = 0 ; i < rv.infiles.size() ; ++i )
    {
      struct stat buf;
      if (stat(rv.infiles[i].c_str(), &buf) == -1) 
	{
	  cerr << "Error: input file "
	       << rv.infiles[i]
	       << " does not exist\n";
	  exit(1);
	}
    }

  return rv;
}

void read_data(putCNVs & raw_div,
	       putCNVs & raw_par,
	       map<string,putCNVs > & raw_ul,
	       const char * filename,
	       const int8_t & min_mqual,
	       const int16_t & max_mm,
	       const int16_t & max_gap)
{
  gzFile lin = gzopen(filename,"r");
  if(lin == NULL)
    {
      cerr << "Error: could not open "
	   << filename
	   << " for reading on line "
	   << __LINE__ << " of "
	   << __FILE__ << '\n';
      exit(1);
    }
  char type[4];
  type[3]='\0';
  do
    {
      auto name = gzreadCstr(lin);
      if(gzeof(lin))break;
      auto chrom = gzreadCstr(lin);
      if( chrom.second <= 0 )
	{
	  cerr << "Error: gzread error on line " << __LINE__
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      auto chrom2 = gzreadCstr(lin);
      if( chrom2.second <= 0 )
	{
	  cerr << "Error: gzread error on line " << __LINE__
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      if( gzread(lin,&type[0],3*sizeof(char)) <= 0 )
	{
	  cerr << "Error: gzread error on line " << __LINE__
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      alnInfo read1(lin),read2(lin);
      if( read1.mapq >= min_mqual && read2.mapq >= min_mqual &&
	  read1.mm <= max_mm && read1.ngap <= max_gap &&
	  read2.mm <= max_mm && read2.ngap <= max_gap )
	{
	  if(string(type) == "DIV")
	    {
	      if ( unique_positions(raw_div[chrom.first],
				    (read1.strand==0) ? read2.start : read1.start,
				    (read1.strand==0) ? read1.start : read2.start) )
		{
		  assert( (read1.strand==0) ? (read2.strand == 1) : (read1.strand == 1) );
		  raw_div[chrom.first].push_back( linkeddata( (read1.strand==0) ? read2.start : read1.start,
							      (read1.strand==0) ? read2.stop : read1.stop,
							      (read1.strand==0) ? read1.start : read2.start,
							      (read1.strand==0) ? read1.stop : read2.stop,
							      name.first,1,0 ) );
		}
	    }
	  else if (string(type) == "PAR")
	    {
	      if ( unique_positions(raw_par[chrom.first],
				    (read1.start<read2.start) ? read1.start : read2.start,
				    (read1.start<read2.start) ? read2.start : read1.start) )
		{
		  raw_par[chrom.first].push_back( linkeddata( (read1.start<read2.start) ? read1.start : read2.start,
							      (read1.start<read2.start) ? read1.stop : read2.stop,
							      (read1.start<read2.start) ? read2.start : read1.start,
							      (read1.start<read2.start) ? read2.stop : read1.stop,
							      name.first,
							      (read1.start<read2.start) ? read1.strand : read2.strand,
							      (read1.start<read2.start) ? read2.strand : read1.strand ));
		}
	    }
	  else if (string(type) == "UNL")
	    {
	      assert(chrom.first != chrom2.first);
	      if( chrom > chrom2 )
		{
		  swap(chrom,chrom2);
		  swap(read1,read2);
		}
	      if ( unique_positions(raw_ul[chrom.first][chrom2.first],read1.start,read2.start) )
		{
		  raw_ul[chrom.first][chrom2.first].push_back( linkeddata(read1.start,read1.stop,
									  read2.start,read2.stop,
									  name.first,
									  read1.strand,read2.strand) );
		}
	    }
#ifndef NDEBUG
	  else
	    {
	      abort();
	    }
#endif
	}
    } while(!gzeof(lin));
  gzclose(lin);
}

bool unique_positions(const lvector & data,
		      const unsigned & start,
		      const unsigned & start2)
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



bool pair_should_cluster( lvector::const_iterator & pair,
			  vector<lvector::const_iterator> & cluster,
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

cluster_container cluster_linked( const lvector & raw,
				  const unsigned & mdist )
{
  using citr = lvector::const_iterator;
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

void write_clusters( gzFile gzout,
			   const string & chrom1,
			   const string & chrom2,
			   const cluster_container & clusters,
			   unsigned * eventid )
{
  for(unsigned i=0;i<clusters.size();++i)
    {
      //get the boundaries of each event
      unsigned min1=numeric_limits<unsigned>::max(),max1=0,min2=numeric_limits<unsigned>::max(),max2=0;
      string readnames;
      for(unsigned j=0;j<clusters[i].size();++j)
	{
	  //The +1 here convert genomic positions to a [1,L] coordinate system
	  min1 = min(min1,clusters[i][j]->a+1);
	  max1 = max(max1,clusters[i][j]->aS+1);
	  min2 = min(min2,clusters[i][j]->b+1);
	  max2 = max(max2,clusters[i][j]->bS+1);
	  ostringstream t;
	  //The +1 here convert genomic positions to a [1,L] coordinate system
	  t << ';' 
	    << clusters[i][j]->a+1 << ',' 
	    << clusters[i][j]->aS+1 << ','
	    << clusters[i][j]->strand1 << ','
	    << clusters[i][j]->b+1 << ',' 
	    << clusters[i][j]->bS+1 << ','
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
      ostringstream o;
      o << *eventid << '\t'
	<< chrom1 << '\t'
	<< clusters[i].size() << '\t'
	<< clusters[i][0]->strand1 << '\t'
	<< min1 << '\t' << max1 << '\t'
	<< chrom2 << '\t'
	<< clusters[i][0]->strand2 << '\t'
	<< min2 << '\t' << max2 << '\t'
	<< readnames << '\n';
      if(!gzwrite(gzout,o.str().c_str(),o.str().size()))
	{
	  cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      ++(*eventid);
    } 
}

void write_clusters_bedpe( gzFile gzout,
			   const string & sampleID,
			   const string & eventtype,
			   const string & chrom1,
			   const string & chrom2,
			   const cluster_container & clusters,
			   unsigned * eventid )
{
  for(unsigned i=0;i<clusters.size();++i)
    {
      //get the boundaries of each event
      unsigned min1=numeric_limits<unsigned>::max(),max1=0,min2=numeric_limits<unsigned>::max(),max2=0;
      string readnames;
      for(unsigned j=0;j<clusters[i].size();++j)
	{
	  //The +1 here convert genomic positions to a [1,L] coordinate system
	  min1 = min(min1,clusters[i][j]->a+1);
	  max1 = max(max1,clusters[i][j]->aS+1);
	  min2 = min(min2,clusters[i][j]->b+1);
	  max2 = max(max2,clusters[i][j]->bS+1);
	  ostringstream t;
	  //The +1 here convert genomic positions to a [1,L] coordinate system
	  t << ';' 
	    << clusters[i][j]->a+1 << ',' 
	    << clusters[i][j]->aS+1 << ','
	    << clusters[i][j]->strand1 << ','
	    << clusters[i][j]->b+1 << ',' 
	    << clusters[i][j]->bS+1 << ','
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
      ostringstream o;
      o << chrom1 << '\t'
	<< (min1-1) << '\t'                                          //b/c start1 is zero-based
	<< max1 << '\t'                                              //b/c start2 is one-based
	<< chrom2 << '\t' 
	<< (min2-1) << '\t'
	<< max2 << '\t'
	//The above are the minimal fields
	<< sampleID << "_" << eventtype << "_event" << *eventid << '\t'                               //"name"
	<< log10(clusters[i].size()) << '\t'                         //The score = log10(coverage)
	<< ( (clusters[i][0]->strand1 == 0 ) ? '+' : '-' ) << '\t'   //strand1
	<< ( (clusters[i][0]->strand2 == 0 ) ? '+' : '-' ) << '\t'   //strand2
	<< readnames <<'\n';                                         //The reads are the optional column
      if(!gzwrite(gzout,o.str().c_str(),o.str().size()))
	{
	  cerr << "Error: gzwrite error encountered at line " << __LINE__ 
	       << " of " << __FILE__ << '\n';
	  exit(1);
	}
      ++(*eventid);
    } 
}
