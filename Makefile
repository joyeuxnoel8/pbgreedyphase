all: partitionByPhasedSNVs readToSNVList

SEQAN=seqan/include
BOOSTLIB=boost/stage/lib
BOOST=boost
BLASR=blasr/common
VCFLIB=vcflib
HTSLIB=$(VCFLIB)/tabixpp/htslib
CPPOPTS=  -g
LZMA=liblzma/dist/lib
CPP=g++ -std=c++14
LIBBZ2=bzip2-1.0.6

vcflib:
	cd vcflib && make -j 8

partitionByPhasedSNVs: PartitionByPhasedSNVs.cpp FastaIndex.h
	$(CPP) -g -static $(CPPOPTS) $< \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -I $(VCFLIB)/include -I $(HTSLIB) \
     -L $(BOOSTLIB) -l boost_program_options \
     -L $(VCFLIB)/lib -lvcflib  -L $(HTSLIB) -l hts -lpthread -lz -L$(LZMA) -llzma -L$(LIBBZ2) -lbz2 -lpthread \
     -o $@ 


readToSNVList: ReadToSNVList.cpp
	$(CPP) -static $(CPPOPTS) $^ \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -L $(BOOSTLIB) -l boost_program_options \
     -o $@ 

