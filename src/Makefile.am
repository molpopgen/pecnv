bin_PROGRAMS=pecnv 

pecnv_SOURCES=pecnv.cc process_readmappings.hpp process_readmappings.cc teclust.cc teclust.hpp common.cc teclust_objects.hpp teclust_objects.cc teclust_phrapify.hpp teclust_phrapify.cc teclust_parseargs.hpp teclust_parseargs.cc teclust_scan_bamfile.hpp teclust_scan_bamfile.cc intermediateIO.hpp intermediateIO.cc cluster_cnv.hpp cluster_cnv2.cc mdist.hpp bwa_mapdistance.cc file_common.hpp file_common.cc mkgenome.hpp mkgenome.cc

AM_CXXFLAGS=
if HAVE_HTSLIB
AM_CXXFLAGS+=-DHAVE_HTSLIB
endif
