CXX=c++
CXXFLAGS=-O2 -Wall -W -ansi -pedantic -I.
SEQ=-lsequence
BIOSTREAMS=-lboost_iostreams
#Modify LDFLAGS if your system has libraries in non-standard locations
targets=fastq_to_table bwa_bam_to_mapfiles bwa_mapdistance rename_reference cluster_cnv teclust umm_te_finder filter_edit get_pair_ids_bwa

all: fastq_to_table.o bwa_bam_to_mapfiles2.o bwa_mapdistance.o bwa_util.o isbinary.o rename_reference.o cluster_cnv2.o string_unsigned_lookup.o file_util.o umm_te_finder.o teclust.o filter_edit.o get_pair_ids_bwa.o
	$(CXX) $(CXXFLAGS) -o rename_reference rename_reference.o isbinary.o $(LDFLAGS) $(SEQ) $(BIOSTREAMS) 
	$(CXX) $(CXXFLAGS) -o fastq_to_table fastq_to_table.o isbinary.o $(LDFLAGS) $(SEQ) $(BIOSTREAMS) 
	$(CXX) $(CXXFLAGS) -o bwa_bam_to_mapfiles bwa_bam_to_mapfiles2.o bwa_util.o isbinary.o $(LDFLAGS) $(SEQ) $(BIOSTREAMS) 
	$(CXX) $(CXXFLAGS) -o bwa_mapdistance bwa_mapdistance.o bwa_util.o isbinary.o $(LDFLAGS) $(SEQ) $(BIOSTREAMS)  
	$(CXX) $(CXXFLAGS) -o cluster_cnv cluster_cnv2.o string_unsigned_lookup.o file_util.o isbinary.o $(LDFLAGS) $(SEQ) $(BIOSTREAMS) 
	$(CXX) $(CXXFLAGS) -o umm_te_finder umm_te_finder.o isbinary.o $(LDFLAGS) $(SEQ) $(BIOSTREAMS) 
	$(CXX) $(CXXFLAGS) -o teclust teclust.o isbinary.o string_unsigned_lookup.o $(LDFLAGS) $(SEQ) $(BIOSTREAMS) 
	$(CXX) $(CXXFLAGS) -o filter_edit filter_edit.o isbinary.o $(LDFLAGS) $(SEQ) $(BIOSTREAMS)
	$(CXX) $(CXXFLAGS) -o get_pair_ids_bwa get_pair_ids_bwa.o isbinary.o $(LDFLAGS) $(SEQ) $(BIOSTREAMS)
dist: 
	mkdir te_alignment_processing
	cp *.cc *.hpp Makefile README.txt te_alignment_processing
	tar czf te_alignment_processing.tar.gz te_alignment_processing
	rm -rf te_alignment_processing
clean:
	rm -f $(targets) *.o

fastq_to_table.o: isbinary.hpp
bwa_bam_to_mapfiles2.o: isbinary.hpp bwa_util.hpp
bwa_mapdistance.o: isbinary.hpp bwa_util.hpp
isbinary.o: isbinary.hpp
bwa_util.o: bwa_util.hpp
rename_reference: isbinary.hpp
cluster_cnv2.o: isbinary.hpp string_unsigned_lookup.o file_util.o
umm_te_finder.o: isbinary.hpp
teclust.o: isbinary.hpp string_unsigned_lookup.hpp
