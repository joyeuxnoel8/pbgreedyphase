all: boost_1_66_0/stage/lib/libboost_program_options.a vcflib/include/Variant.h partitionByPhasedSNVs readToSNVList 

SEQAN=seqan/include
BOOST=boost_1_66_0
BOOSTLIB=$(BOOST)/stage/lib
BLASR=blasr/common
VCFLIB=vcflib
HTSLIB=$(VCFLIB)/tabixpp/htslib
CPPOPTS=  -g
LZMA=liblzma/build/lib
CPP=g++ -std=c++14
LIBBZ2=bzip2-1.0.6

vcflib/include/Variant.h:
	cd vcflib && make -j 8

boost_1_66_0.tar.gz:
	wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz
	tar xvf boost_1_66_0.tar.gz

boost_1_66_0/stage/lib/libboost_program_options.a: boost_1_66_0.tar.gz
	cd boost_1_66_0 && ./bootstrap.sh && ./b2 --prefix=$PWD/build -j 4

partitionByPhasedSNVs: PartitionByPhasedSNVs.cpp FastaIndex.h boost_1_66_0/stage/lib/libboost_program_options.a
	$(CPP) -g -static $(CPPOPTS) $< \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -I $(VCFLIB)/include -I $(HTSLIB) \
     -L $(BOOSTLIB) -l boost_program_options \
     -L $(VCFLIB)/lib -lvcflib  -L $(HTSLIB) -l hts -lpthread -lz -L$(LZMA) -llzma -L$(LIBBZ2) -lbz2 -lpthread \
     -o $@ 


readToSNVList: ReadToSNVList.cpp PartitionTools.h FastaIndex.h SamUtils.h GenotypedRead.h SNVDB.h 
	$(CPP) -static $(CPPOPTS) $< \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -L $(BOOSTLIB) -l boost_program_options \
     -o $@ 

