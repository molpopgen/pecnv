#include <teclust_parseargs.hpp>
#include <boost/program_options.hpp>
#include <sys/stat.h>
using namespace std;
using namespace boost::program_options;

teclust_params teclust_parseargs(const int argc, char ** argv)
{
  teclust_params rv;
  options_description desc("Cluster reads into putative transposable element calls.\nUsage: teclust -h to see help");
  desc.add_options()
    ("help,h", "Produce help message")
    ("bamfile,b",value<string>(&rv.bamfile),"BAM file name (optional)")
    ("tepos,t",value<string>(&rv.reference_datafile),"File containing positions of TEs in reference genome (optional)")
    ("outfile,o",value<string>(&rv.outfile),"Output file name for clusters (required)")
    ("umu,u",value<string>(&rv.umufile),"The um_u output file for the sample generated by process_readmappings (required)")
    ("umm,m",value<string>(&rv.ummfile),"The um_u output file for the sample generated by process_readmappings (required)")
    //("umuref,r",value<string>(&rv.teclust_ref),"The output of this program run on a reference genome, if available (optional)")
    ("isize,i",value<int32_t>(&rv.INSERTSIZE),"Upper limit on insert size distribution, e.g. from bwa_mapdistance. (required)")
    ("mdist,M",value<int32_t>(&rv.MDIST)->default_value(1000),"Max. distance for joining up left and right ends of putative TEs (required)")
    ("phrapdir,p",value<string>(&rv.phrapdir),"Name of a directory to put input files for de novo assembly of putatitve TE insertions using phrap. If the directory does not exist, it will be created. (optional)")
    ("minreads,r",value<int32_t>(&rv.MINREADS)->default_value(3),"Min. number of reads in a cluster for writing input files for phrap. (optional)")
    ("closestTE,c",value<int>(&rv.CLOSEST)->default_value(-1),"For phrap output, only consider events >= c bp away from closest TE in the reference. (optional)")
    ("ummHitTE","When processing the um_u/um_m files, only consider reads where the M read hits a known TE.  This makes --tepos/-t a required option. (optional)")
    ("allEvents,a","For phrap output: write files for all events. Default is only to write files for putative novel insertions");
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

  if(vm.count("allEvents"))
    {
      rv.novelOnly=false;
    }

  if(vm.count("ummHitTE"))
    {
      rv.greedy = false;
      if(! vm.count("tepos"))
	{
	  cerr << "Error: --tepos/-t is required if --ummHitTe is specified.\n"
	       << desc << '\n';
	  exit(0);
	}
    }

  //Make sure that values are sane
  if( rv.INSERTSIZE <= 0 ) 
    {
      cerr << "Error: value passed to --isize/-i must be > 0\n";
      exit(0);
    }
  if( rv.MDIST <= 0 ) 
    {
      cerr << "Error: value passed to --mdist/-M must be > 0\n";
      exit(0);
    }
  if( rv.MINREADS <= 0 ) 
    {
      cerr << "Error: value passed to --minreads/-r must be > 0\n";
      exit(0);
    }

  //Check that specified input files exist
  struct stat buf;
  if( vm.count("bamfile") )
    {
      if (stat(rv.bamfile.c_str(), &buf) == -1) 
	{
	  cerr << "Error: "
	       << rv.bamfile
	       << " does not exist\n";
	}
    }
  if ( vm.count("tepos") )
    {
      if (stat(rv.reference_datafile.c_str(), &buf) == -1) 
	{
	  cerr << "Error: "
	       << rv.reference_datafile
	       << " does not exist\n";
	}
    }
  if (stat(rv.umufile.c_str(), &buf) == -1) 
    {
      cerr << "Error: "
	   << rv.umufile
	   << " does not exist\n";
    }
  if (stat(rv.ummfile.c_str(), &buf) == -1) 
    {
      cerr << "Error: "
	   << rv.ummfile
	   << " does not exist\n";
    }
  return rv;
}
