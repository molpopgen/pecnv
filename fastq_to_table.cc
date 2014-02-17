//reformat 1 fastq file from stdin

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <cctype>
//this project
#include <isbinary.hpp>
//libsequence
#include <Sequence/Fasta.hpp>
//boost
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>


using namespace std;
using namespace boost::iostreams;

template<typename streamtype>
void process_files(streamtype & in,
		   const unsigned line_id,
		   const unsigned lane_id,
		   const unsigned read_id,
		   //const char * ofn,
		   const char * ofn2,
		   const bool & revcom)
{
  string name,seq,name2,qual,qual2;

  unsigned pair_id=0;

  //ofstream os(ofn,ios_base::out | ios_base::binary);
  //filtering_ostream o;
  //o.push(gzip_compressor());
  //o.push(os);
  //o.push(file_sink(ofn),ios_base::out | ios_base::binary);
  //ofstream os2(ofn2,ios_base::out | ios_base::binary);
  filtering_ostream o2;
  o2.push(gzip_compressor());
  //o2.push(os2);
  o2.push(file_sink(ofn2),ios_base::out | ios_base::binary);
  while(! in.eof() )
    {
      getline(in,name);
      getline(in,seq);
      getline(in,name2);
      getline(in,qual);
      in >> ws;
      string::iterator i,j=name.begin();
      i = find_if(name.begin(),name.end(),::isspace);
      while(i != name.end())
	{
	  name.erase(i);
	  j=i+1;
	  i=find_if(j,name.end(),::isspace);
	}

      i = find_if(name2.begin(),name2.end(),::isspace);
      while(i != name2.end())
	{
	  name2.erase(i);
	  j=i+1;
	  i=find_if(j,name2.end(),::isspace);
	}


      if( revcom )
	{
	  reverse(qual.begin(),qual.end());
	  Sequence::Fasta temp("temp",seq);
	  temp.Revcom();
	  seq = temp.second;
	}
      
      o2 <<'@'<< line_id << ':' << lane_id << ':' << pair_id << ':' << read_id <<'\n'
	 << seq << '\n'
	 << '+' << line_id << ':' << lane_id << ':' << pair_id << ':' << read_id <<'\n'
	 << qual << endl;
      
      /*
      o << line_id << '\t' << lane_id << '\t' << pair_id++ << '\t'
	<< read_id << '\t' << name << '\t' << seq << '\t' << qual << endl;
      */
    }
}

int main(int argc, char ** argv)
{
  int argn=1;
  if( argc < 7 )
    {
      cerr << "usage: " << argv[0]
	   << " line_id lane_id read_id fastqfile table_file_name new_fastq_file_name\n";
      exit(10);
    }
  const unsigned line_id = atoi(argv[argn++]);
  const unsigned lane_id = atoi(argv[argn++]);
  const unsigned read_id = atoi(argv[argn++]);
  const char * fastqfile = argv[argn++];
  //const char * ofn =  argv[argn++];
  const char * ofn2 = argv[argn++];
  bool revcom = false;
  if( argc == 8 )
    {
      revcom = atoi(argv[argn++]);
    }
  if ( isbz2(fastqfile) )
    {
      filtering_istream in;
      in.push(bzip2_decompressor());
      in.push(file_source(fastqfile,ios_base::in|ios_base::binary));
      //process_files(in,line_id,lane_id,read_id,ofn,ofn2,revcom);
      process_files(in,line_id,lane_id,read_id,ofn2,revcom);
    }
  else if ( isbinary(fastqfile) )
    {
      filtering_istream in;
      in.push(gzip_decompressor());
      in.push(file_source(fastqfile,ios_base::in|ios_base::binary));
      process_files(in,line_id,lane_id,read_id,ofn2,revcom);
    }
  else
    {
      ifstream in(fastqfile);
      process_files(in,line_id,lane_id,read_id,ofn2,revcom);
    }
}



