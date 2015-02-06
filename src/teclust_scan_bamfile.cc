#include <teclust_scan_bamfile.hpp>
#include <Sequence/bamreader.hpp>
#include <Sequence/bamrecord.hpp>
#include <Sequence/samfunctions.hpp>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <common.hpp>
#include <sys/stat.h>

using namespace std;
using namespace Sequence;

using puu = pair<int32_t,int8_t>;

using refIDlookup = unordered_map<int32_t,string>;

//DEFINITION OF FUNCTIONS
refIDlookup make_lookup(const bamreader & reader);

void scan_bamfile(const teclust_params & p,
		  const refTEcont & refTEs,
		  unordered_set<string> * readPairs,
		  map<string,vector< puu > > * data,
		  hts_idx_t * idx)
{
  if( refTEs.empty() || p.bamfile.empty() ) return; 
  struct stat buf;
  if (stat(p.bamfile.c_str(), &buf) == -1) 
    {
      cerr << "Error: "
	   << p.bamfile
	   << " does not exist\n";
    }
  bamreader reader(p.bamfile.c_str());
  if(! reader )
    {
      cerr << "Error: " 
	   << p.bamfile
	   << " could not be opened for reading.\n";
      exit(0);
    }

  auto lookup = make_lookup(reader);

  auto firstREC = reader.tell();
  //unordered_set<string> RPlocal;
  unordered_map< string, pair<int32_t,int32_t> > RPlocal;
  while(! reader.eof() && !reader.error() )
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break;
      samflag f(b.flag());
      if( !f.query_unmapped && !f.mate_unmapped )
	//then both reads are mapped 
	{
	  auto n = editRname(b.read_name());
	  if( readPairs->find(n) == readPairs->end())
	    {
	      auto itr = lookup.find(b.refid());
	      
	      if(itr == lookup.end())
		{
		  cerr << "Error: reference ID " << b.refid()
		       << " not found in BAM file header. Line "
		       << __LINE__ 
		       << " of " << __FILE__ << '\n';
		  exit(1);
		}
	      
	      //Now, does the read overlap a known TE?
	      int32_t start = b.pos(),stop=b.pos() + alignment_length(b) - 1;
	      auto CHROM = refTEs.find(itr->second);
	      if( CHROM != refTEs.end() )
		{
		  bool hitsTE = find_if( CHROM->second.cbegin(),
					 CHROM->second.cend(),
					 [&](const teinfo & __t) {
					   bool A = (start >= __t.start() && start <= __t.stop());
					   bool B = (stop >= __t.start() && stop <= __t.stop());
					   return A||B;
					 }) != CHROM->second.cend();
		  if( hitsTE )
		    {
		      /*We can do a check here:
			If mate is mapped to same chromo & hits a TE,
			we can skip storing it.
			We don't have access to it's mates start and stop,
			but we can check the start.
		      */
		      bool OK = true;
		      if(b.refid() == b.next_refid())
			{
			  int32_t mstart = b.next_pos();
			  OK = find_if( CHROM->second.cbegin(),
					CHROM->second.cend(),
					[&](const teinfo & __t) {
					  return (mstart >= __t.start() && mstart <= __t.stop());
					}) == CHROM->second.cend();
			}
		      if(OK)
			{
			  //readPairs->insert( n );
			  //RPlocal.insert(n);
			  RPlocal[n] = make_pair( b.next_refid(), b.next_pos() );
			}
		    }
		}
	    }
	}
    }

  //SECOND PASS -- new version based on seeking within the bgzf file
  cerr << "Second pass starting.  There are " << RPlocal.size() << " names to process\n";
  unsigned FOUND = 0;
  for( auto itr = RPlocal.cbegin(); itr != RPlocal.cend() ; ++itr )
    {
      /*Step 1: seek to where the record should be.
	We are seeking to the position where this read should be to this position + 1
	htslib doesn't seem to like to seek to an i,i range.
       */
      hts_itr_t *iter = bam_itr_queryi(idx,itr->second.first,itr->second.second,itr->second.second+1);
      if(iter->off != NULL) //If there are reads in the range
	{
	  reader.seek(iter->off->u,SEEK_SET); //seek to position
	  bool found = false;
	  while(!found && reader.tell() < iter->off->v )
	    {
	      bamrecord b = reader.next_record();
	      samflag f(b.flag());
	      if( !f.query_unmapped && !f.mate_unmapped )
		{
		  auto itr2 = lookup.find(b.refid());
		  string n = editRname(b.read_name());
		  if(n == itr->first)
		    {
		      found = true;
		      ++FOUND;
		      if( FOUND > 0 && FOUND % 100 == 0. ) cerr << FOUND << '\n';
		      int32_t start = b.pos(),stop=b.pos() + alignment_length(b) - 1;
		      auto CHROM = refTEs.find(itr2->second);
		      if( CHROM != refTEs.end() )
			{
			  bool hitsTE = find_if( CHROM->second.cbegin(),
						 CHROM->second.cend(),
						 [&](const teinfo & __t) {
						   bool A = (start >= __t.start() && start <= __t.stop());
						   bool B = (stop >= __t.start() && stop <= __t.stop());
						   return A||B;
						 }) != CHROM->second.cend();
			  if(!hitsTE)
			    {
			      auto DCHROM = data->find(itr2->second);
			      if(DCHROM == data->end())
				{
				  data->insert( make_pair(itr2->second,
							  vector<pair<int32_t,int8_t> >(1,make_pair(start,f.qstrand) ) ) );
				}
			      else
				{
				  DCHROM->second.emplace_back(make_pair(start,f.qstrand));
				}
			    }
			}
		    }
		}
	    }
	}
      //cleanup
      hts_itr_destroy(iter);
    }
  //Second pass -- OLD VERSION
  // reader.seek( firstREC, SEEK_SET );
  // while(! reader.eof() && !reader.error() )
  //   {
  //     bamrecord b = reader.next_record();
  //     if(b.empty()) break;
  //     samflag f(b.flag());
  //     if( !f.query_unmapped && !f.mate_unmapped) 
  // 	//then both reads are mapped 
  // 	{
  // 	  auto itr = lookup.find(b.refid());
  // 	  if(itr == lookup.end())
  // 	    {
  // 	      cerr << "Error: reference ID " << b.refid()
  // 		   << " not found in BAM file header. Line "
  // 		   << __LINE__ 
  // 		   << " of " << __FILE__ << '\n';
  // 	      exit(1);
  // 	    }
  // 	  auto n = editRname(b.read_name());
  // 	  /*
  // 	    Note: Julie's script does not check that these reads map uniquely.
  // 	  */
  // 	  //if( readPairs->find(n) != readPairs->end() )
  // 	  if( RPlocal.find(n) != RPlocal.end() && readPairs->find(n) == readPairs->end() )
  // 	    {
  // 	      int32_t start = b.pos(),stop= b.pos() + alignment_length(b) - 1;
  // 	      auto CHROM = refTEs.find(itr->second);
  // 	      if( CHROM != refTEs.end() )
  // 		{
  // 		  bool hitsTE = find_if( CHROM->second.cbegin(),
  // 					 CHROM->second.cend(),
  // 					 [&](const teinfo & __t) {
  // 					   bool A = (start >= __t.start() && start <= __t.stop());
  // 					   bool B = (stop >= __t.start() && stop <= __t.stop());
  // 					   return A||B;
  // 					 }) != CHROM->second.cend();
  // 		  if(!hitsTE)
  // 		    {
  // 		      auto DCHROM = data->find(itr->second);
  // 		      if(DCHROM == data->end())
  // 			{
  // 			  data->insert( make_pair(itr->second,
  // 						  vector<pair<int32_t,int8_t> >(1,make_pair(start,f.qstrand) ) ) );
  // 			}
  // 		      else
  // 			{
  // 			  DCHROM->second.emplace_back(make_pair(start,f.qstrand));
  // 			}
  // 		    }
  // 		}
  // 	    }
  // 	}
  //   }
}

refIDlookup
make_lookup(const bamreader & reader)
{
  refIDlookup rv;
  unsigned I=0;
  for_each(reader.ref_cbegin(),reader.ref_cend(),
	   [&](const pair<string,int32_t> & __p)
	   {
	     rv[I++]=__p.first;
	   });
  return rv;
}
