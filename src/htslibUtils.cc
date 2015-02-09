#include <htslibUtils.hpp>
#include <htslib/bgzf.h>
#include <sys/stat.h>

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
