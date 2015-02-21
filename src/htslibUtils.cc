#include <htslibUtils.hpp>
#include <htslib/bgzf.h>
#include <sys/stat.h>
#include <limits>
#include <map>

using namespace std;

hts_idx_t * get_index(const std::string & bamfilename)
{
  string indexfilename = bamfilename + ".bai";
  struct stat buf;
  if (stat(indexfilename.c_str(), &buf) == -1) return nullptr;
  BGZF * bam = bgzf_open(bamfilename.c_str(),"rb");
  hts_idx_t * idx = bam_index_load(indexfilename.c_str());
  bgzf_close(bam);
  return idx;
}

//use htslib to find out where chromosomes begin and end in our bam file
//Christ, they really need to document that fucking library...
vector<pair<string,pair<uint64_t,uint64_t> >> read_index(const std::string & bamfilename)
{
  string indexfilename = bamfilename + ".bai";

  vector<pair<string,pair<uint64_t,uint64_t> > > rv;
  //open fam file, read header
  BGZF * bam = bgzf_open(bamfilename.c_str(),"rb");
  bam_hdr_t * hdr = bam_hdr_read(bam);


  //read in the index
  hts_idx_t * idx = bam_index_load(indexfilename.c_str());
  if( idx == nullptr )
    {
      bam_hdr_destroy(hdr);
      bgzf_close(bam);
      return rv;
    }
  for( int32_t i = 0 ; i < hdr->n_targets ; ++i )
    {
      hts_itr_t *iter = bam_itr_queryi(idx,i,0,hdr->target_len[i]);
      if(iter->off != NULL )
	{
	  rv.push_back( make_pair( string(hdr->target_name[i]),
				   make_pair(iter->off->u,iter->off->v) ));
	}
      hts_itr_destroy(iter);
    }
  //cleanup
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  bgzf_close(bam);
  return rv;
} 

bamrange::bamrange( const std::vector<refrange> & ranges ) : beg( numeric_limits<uint64_t>::max() ),
							     end( numeric_limits<uint64_t>::min() ),
							     refid1(-1),
							     refid2(-1),
							     start(0),stop(0),
							     refname1(std::string()),
							     refname2(std::string())
{
  for( const refrange & r : ranges )
    {
      if( r.beg < beg ) {
	beg = r.beg;
	refid1 = r.refid;
	start = r.start;
	refname1 = r.refname;
      }
      if( r.end > end ) {
	end = r.end;
	refid2 = r.refid;
	stop = r.stop;
	refname2 = r.refname;
      }
    }
}

vector<bamrange> split_genome(const hts_idx_t * idx,
			      const std::string & bamfilename,
			      const unsigned & NTHREADS)
{
  if(idx == nullptr) return vector<bamrange>();

  BGZF * bam = bgzf_open(bamfilename.c_str(),"rb");
  bam_hdr_t * hdr = bam_hdr_read(bam);

  unsigned genome_len = 0;
  vector<pair<pair<string,int32_t>,pair<unsigned,unsigned> > >lengths;
  for( int32_t i = 0 ; i < hdr->n_targets ; ++i )
    {    
      hts_itr_t *iter = bam_itr_queryi(idx,i,0,hdr->target_len[i]);
      if(iter->off != NULL )
	{
	  lengths.push_back(make_pair(make_pair(hdr->target_name[i],i),
				      make_pair(1,hdr->target_len[i])));
	  genome_len += hdr->target_len[i];
	}      
      hts_itr_destroy(iter);
    }
  unsigned chunksize = genome_len/(NTHREADS),chunkmainder = unsigned(genome_len%NTHREADS);
  
  unsigned CC=chunksize+chunkmainder;
  unsigned j = 0;
  
  //Need to store these results in some container...
  map<unsigned,vector< refrange > > chunks;
  for(unsigned i = 0 ; j<lengths.size()&&i < NTHREADS;)
    {
      for( ; CC>0 && j < lengths.size() ;  )
	{
	  auto CLEN = lengths[j].second.second - lengths[j].second.first + 1;
	  if( CLEN <= CC )
	    {      
	      hts_itr_t *iter = bam_itr_queryi(idx,lengths[j].first.second,lengths[j].second.first-1,lengths[j].second.second);
	      chunks[i].emplace_back( refrange(iter->off->u, iter->off->v,
					       lengths[j].first.second,lengths[j].second.first-1,lengths[j].second.second,lengths[j].first.first) );
	      hts_itr_destroy(iter);
	      CC -= CLEN;
	      ++j;
	    }
	  else
	    {
	      hts_itr_t *iter = bam_itr_queryi(idx,lengths[j].first.second,lengths[j].second.first-1,lengths[j].second.first-1+CC);
	      chunks[i].emplace_back( refrange(iter->off->u, iter->off->v,
					       lengths[j].first.second,
					       lengths[j].second.first-1,
					       lengths[j].second.first-1+CC,
					       lengths[j].first.first) );
	      hts_itr_destroy(iter);
	      
	      lengths[j].second.first += CC;
	      ++i;
	      CC=chunksize;
	    }
	}
    }
  vector< bamrange > ranges(NTHREADS);
  for ( const auto & __p : chunks )
    {
      ranges[__p.first]=bamrange(__p.second);
    }
  bam_hdr_destroy(hdr);
  bgzf_close(bam);
  return ranges;
}
