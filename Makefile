all: boost_1_66_0/stage/lib/libboost_program_options.a vcflib_build_flag partitionByPhasedSNVs readToSNVList 

SEQAN=seqan/include
BOOST=boost_1_66_0
BOOSTLIB=$(BOOST)/stage/lib
BLASR=blasr/common
VCFLIB=vcflib
HTSLIB=$(VCFLIB)/tabixpp/htslib
CPPOPTS=  -g
LZMA=lzma/build/lib
CPP=g++ -std=c++14
LIBBZ2=bzip2-1.0.6

VCFLIB_INCLUDES := "-I$(abspath vcflib/tabixpp/htslib) -I$(abspath vcflib/include) -L. -L$(abspath vcflib/tabixpp/htslib) -I$(abspath $(LIBBZ2)) -L$(abspath $(LIBBZ2))"
VCFLIB_LDFLAGS := "-Llib -lvcflib -lhts -lpthread -lz -lm -llzma -lbz2"
#VCFLIB_FLAGS = "-I$(abspath $(LIBBZ2))"

$(LIBBZ2)/libbz2.a:
	cd $(LIBBZ2) && make

$(LZMA)/liblzma.a:
	cd liblzma && ./configure --prefix=$(PWD)/lzma/build && make && make install

vcflib_build_flag: $(LIBBZ2)/libbz2.a
	make -C vcflib -j 8 INCLUDES=$(VCFLIB_INCLUDES) LDFLAGS=$(VCFLIB_LDFLAGS) CFLAGS=$(VCFLIB_INCLUDES)
	touch $@

boost_1_66_0/bootstrap.sh:
	wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz
	tar xvf boost_1_66_0.tar.gz
	touch $@

boost_1_66_0/stage/lib/libboost_program_options.a: boost_1_66_0/bootstrap.sh
	cd boost_1_66_0 && ./bootstrap.sh --without-libraries=python && ./b2 --prefix=$PWD/build -j 4

partitionByPhasedSNVs: PartitionByPhasedSNVs.cpp FastaIndex.h boost_1_66_0/stage/lib/libboost_program_options.a $(LIBBZ2)/libbz2.a $(LZMA)/liblzma.a
	$(CPP) -g $(CPPOPTS) -static $^ \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -I $(VCFLIB)/include -I $(HTSLIB) \
     -L $(BOOSTLIB) -l boost_program_options \
     -L $(VCFLIB)/lib -lvcflib  -L $(HTSLIB) -l hts -lpthread -lz -L$(LIBBZ2) -lbz2 -lpthread -L$(LZMA) -l lzma \
     -o $@ 


readToSNVList: ReadToSNVList.cpp PartitionTools.h FastaIndex.h SamUtils.h GenotypedRead.h SNVDB.h 
	$(CPP) $(CPPOPTS) $< \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -L $(BOOSTLIB) -l boost_program_options \
     -o $@ 

.PHONY: clean
clean:
	make -C $(LIBBZ2) -f Makefile-libbz2_so clean
	rm $(LIBBZ2)/libbz2.so
	make -C vcflib clean
	rm -rf boost_1_66_0.tar.gz boost_1_66_0 partitionByPhasedSNVs readToSNVList
