#bin_PROGRAMS=fastq_to_table bwa_bam_to_mapfiles bwa_mapdistance rename_reference cluster_cnv teclust umm_te_finder filter_edit process_readmappings 
bin_PROGRAMS=cluster_cnv teclust process_readmappings bwa_mapdistance

#get_pair_ids_bwa

#fastq_to_table_SOURCES=fastq_to_table.cc isbinary.cc
#bwa_bam_to_mapfiles_SOURCES=bwa_bam_to_mapfiles2.cc bwa_util.cc
bwa_mapdistance_SOURCES=bwa_mapdistance.cc common.cc common.hpp
#rename_reference_SOURCES=rename_reference.cc isbinary.cc
cluster_cnv_SOURCES=cluster_cnv2.cc intermediateIO.hpp intermediateIO.cc #string_unsigned_lookup.cc file_util.cc file_util.hpp string_unsigned_lookup.hpp
teclust_SOURCES=teclust.cc common.cc teclust_objecst.hpp teclust_objects.cc teclust_phrapify.hpp teclust_phrapify.cc teclust_parseargs.hpp teclust_parseargs.cc teclust_scan_bamfile.hpp teclust_scan_bamfile.cc intermediateIO.hpp intermediateIO.cc #string_unsigned_lookup.cc string_unsigned_lookup.hpp
#umm_te_finder_SOURCES=umm_te_finder.cc 
#filter_edit_SOURCES=filter_edit.cc
process_readmappings_SOURCES=process_readmappings.cc common.cc common.hpp intermediateIO.cc intermediateIO.hpp
#get_pair_ids_bwa_SOURCES=get_pair_ids_bwa.cc isbinary.cc

AM_CXXFLAGS=
if HAVE_HTSLIB
AM_CXXFLAGS+=-DHAVE_HTSLIB
endif
