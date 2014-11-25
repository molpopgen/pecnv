#include <pecnv_version.hpp>
#include <teclust.hpp>
#include <process_readmappings.hpp>
#include <cluster_cnv.hpp>
#include <mdist.hpp>
#include <mkgenome.hpp>
#include <algorithm>
#include <cstring>
#include <iostream>

using namespace std;

void usage(int status);
void citation(int status);
char ** strip_argv(int argc, char ** argv, const char * pattern);

int main( int argc, char ** argv )
{
  if( argc < 2 )
    {
      usage(0);
    }

  //What module will re run?
  if ( strcmp(argv[1],"help") == 0 )
    {
      usage(0);
    }
  if ( strcmp(argv[1],"version") == 0 )
    {
      cout << "pecnv version " << PECNV_VERSION << '\n';
      exit(0);
    }
  else if ( strcmp(argv[1],"process") == 0 )
    {
      auto x = strip_argv(argc,argv,argv[1]);
      process_readmappings_main(x - argv, argv);
    }
  else if( strcmp(argv[1],"cnvclust") == 0 )
    {
      auto x = strip_argv(argc,argv,argv[1]);
      cluster_cnv_main(x - argv, argv);
    }
  else if( strcmp(argv[1],"mdist") == 0 )
    {
      auto x = strip_argv(argc,argv,argv[1]);
      bwa_mapdistance_main(x - argv, argv);
    }
  else if( strcmp(argv[1],"teclust") == 0 )
    {
      auto x = strip_argv(argc,argv,argv[1]);
      teclust_main(x - argv, argv);
    }
  else if( strcmp(argv[1],"mkgenome") == 0 )
    {
      auto x = strip_argv(argc,argv,argv[1]);
      mkgenome_main(x - argv,argv);
    }
  else if( strcmp(argv[1],"citation") == 0 )
    {
      citation(0);
    }
  else
    {
      cerr << "Unkown module: " << argv[1] << '\n';
      usage(0);
    }
}

void usage(int status)
{
  cerr << "pecnv version " << PECNV_VERSION << '\n'
       << "Usage: pecnv module\n"
       << "Available modules are available for:\n"
       << "\tProcessing data in BAM files:\n"
       << "\t\tprocess - collect unusual read pairs from BAM file\n"
       << "\t\tmdist - estimate insert size distribution from BAM file\n"
       << "\tClustering data collected from BAM files:\n"
       << "\t\tcnvclust - perform CNV clustering based on results from process step\n"
       << "\t\tteclust - perform TE clustering based on results from process step\n"
       << "\tManipulating data files\n"
       << "\t\tmkgenome - makes a \"genome file\" for bedtools from a fasta input file\n"
       << "\tProviding info about this program:\n"
       << "\t\tversion - print version info to stdout\n"
       << "\t\tcitation - print citation info to stdout\n";
  exit(status);
}

char ** strip_argv(int argc, char ** argv, const char * pattern)
{
  return remove_if(argv,argv+argc,[&](const char * __xx) { return std::strcmp(__xx,pattern) == 0; });
}

void citation(int status)
{
  cout << "Citation info for pecnv version " << PECNV_VERSION << "\n\n"
       << "The software is available from http://github.com/molpopgen/pecnv\n\n"
       << "The CNV clustering methods are described in:\n\n"
       << "\tRogers, R. L.,J. M. Cridland, L. Shao, T. T. Hu, P. Andolfatto and K. R. Thornton (2014) Landscape of standing variation for tandem duplications in Drosophila yakuba and Drosophila simulans.  Molecular Biology and Evolution 31: 1750-1766 PMID 24710518\n\n"
       << "The TE detection methods are described in:\n\n"
       << "\tCridland, J.M., S.J. MacDonald, A.D. Long, and K.R Thornton (2013) Abundance and Distribution of Transposable Elements in Two Drosophila QTL Mapping Resources  Molecular Biology and Evolution 30: 2311-2327. PMID 23883524"
       << endl;
  exit(status);
}
