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
ZLIB=zlib

$(ZLIB)/build/lib/libz.a:
	cd $(ZLIB) && ./configure --static --prefix=$(PWD)/zlib/build && make -j 8 && make install

vcflib/include/Variant.h:
	cd vcflib && make -j 8

boost_1_66_0/bootstrap.sh:
	wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz
	tar xvf boost_1_66_0.tar.gz
	touch $@

boost_1_66_0/stage/lib/libboost_program_options.a: boost_1_66_0/bootstrap.sh
	cd boost_1_66_0 && ./bootstrap.sh --without-libraries=python && ./b2 --prefix=$PWD/build -j 4

partitionByPhasedSNVs: PartitionByPhasedSNVs.cpp FastaIndex.h boost_1_66_0/stage/lib/libboost_program_options.a $(ZLIB)/build/lib/libz.a
	$(CPP) -g $(CPPOPTS) $< \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -I $(VCFLIB)/include -I $(HTSLIB) \
     -L $(BOOSTLIB) -l boost_program_options \
     -L $(VCFLIB)/lib -lvcflib  -L $(HTSLIB) -l hts -lpthread -L $(ZLIB)/build/lib -lz -L$(LZMA) -llzma -L$(LIBBZ2) -lbz2 -lpthread \
     -o $@ 


readToSNVList: ReadToSNVList.cpp PartitionTools.h FastaIndex.h SamUtils.h GenotypedRead.h SNVDB.h 
	$(CPP) $(CPPOPTS) $< \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -L $(BOOSTLIB) -l boost_program_options \
     -o $@ 

