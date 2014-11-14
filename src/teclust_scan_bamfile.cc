#include <teclust_scan_bamfile.hpp>
#include <Sequence/bamreader.hpp>
#include <Sequence/bamrecord.hpp>
#include <Sequence/samfunctions.hpp>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <common.hpp>
using namespace std;
using namespace Sequence;

using puu = pair<unsigned,unsigned>;

using refIDlookup = unordered_map<int32_t,string>;

//DEFINITION OF FUNCTIONS
refIDlookup make_lookup(const bamreader & reader);

unordered_set<string> scan_bamfile(const params & p,
				   const refTEcont & refTEs,
				   map<string,vector< puu > > * data)
{
  if( refTEs.empty() || p.bamfile.empty() ) return unordered_set<string>();

  bamreader reader(p.bamfile.c_str());
  if(! reader )
    {
      cerr << "Error: " 
	   << p.bamfile
	   << " could not be opened for reading.\n";
      exit(0);
    }

  auto lookup = make_lookup(reader);

  unordered_set<string> TEhitters;
  auto firstREC = reader.tell();
  while(! reader.eof() && !reader.error() )
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break;
      samflag f(b.flag());
      if( !f.query_unmapped && !f.mate_unmapped )
	//then both reads are mapped 
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
	  int32_t start = b.pos(),stop=b.pos() + alignment_length(b);
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
		  //We can do a check here:
		  //If mate is mapped to same chromo & hits a TE,
		  //we can skip storing it
		  bool OK = true;
		  if( b.refid() == b.next_refid())
		    {
		      int32_t mstart = b.next_pos();
		      OK = find_if( CHROM->second.cbegin(),
				    CHROM->second.cend(),
				    [&](const teinfo & __t) {
				      return (mstart >= __t.start() && mstart <= __t.stop());
				}) == CHROM->second.cend();
		    }
		  if(OK)
		    TEhitters.insert( editRname(b.read_name()) );
		}
	    }
	}
    }
  //Second pass
  reader.seek( firstREC, SEEK_SET );
  while(! reader.eof() && !reader.error() )
    {
      bamrecord b = reader.next_record();
      if(b.empty()) break;
      samflag f(b.flag());
      if( !f.query_unmapped && !f.mate_unmapped )
	//then both reads are mapped 
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
	  auto n = editRname(b.read_name());
	  /*
	    Note: Julie's script does not check that these reads map uniquely.
	    Here, we do.
	  */
	  if( TEhitters.find(n) != TEhitters.end() )
	    {
	      //Our check for uniqueness is the existence of an XO flag with a value of 0
	      bamaux ba = b.aux("XO");
	      if(ba.size && ba.value[0] == '1')
		{
		  int32_t start = b.pos(),stop=b.pos() + alignment_length(b);
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
		      if(hitsTE) //no good
			{
			  TEhitters.erase(n);
			}
		      else
			{
			  auto DCHROM = data->find(itr->second);
			  if(DCHROM == data->end())
			    {
			  data->insert( make_pair(itr->second,
						  vector<pair<unsigned,unsigned> >(1,make_pair(start,f.qstrand) ) ) );
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
  return TEhitters;
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
