all: boost_1_66_0/stage/lib/libboost_program_options.a vcflib/lib/libvcflib.a partitionByPhasedSNVs readToSNVList 

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
ZLIB=zlib

$(ZLIB)/build/lib/libz.a:
	cd $(ZLIB) && ./configure --prefix=$(PWD)/zlib/build && make -j 8 && make install

VCFLIB_INCLUDES := "-I$(abspath vcflib/tabixpp/htslib) -I$(abspath vcflib/include) -L. -L$(abspath vcflib/tabixpp/htslib) -I$(abspath $(LIBBZ2)) -L$(abspath $(LIBBZ2))"
VCFLIB_LDFLAGS := "-Llib -lvcflib -lhts -lpthread -lz -lm -L $(LZMA) -llzma -lbz2"


$(LIBBZ2)/libbz2.a:
	cd $(LIBBZ2) && make

$(LZMA)/liblzma.a:
	cd liblzma && ./configure --disable-nls --prefix=$(PWD)/lzma/build && make -j 8 && make install

vcflib/lib/libvcflib.a: $(LIBBZ2)/libbz2.a $(LZMA)/liblzma.a
	make -C vcflib -j 8 INCLUDES=$(VCFLIB_INCLUDES) LDFLAGS=$(VCFLIB_LDFLAGS) CFLAGS=$(VCFLIB_INCLUDES)
	touch $@

boost_1_66_0/bootstrap.sh:
	wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz
	tar xvf boost_1_66_0.tar.gz
	touch $@

boost_1_66_0/stage/lib/libboost_program_options.a: boost_1_66_0/bootstrap.sh
	cd boost_1_66_0 && ./bootstrap.sh --with-libraries=program_options && ./b2 --prefix=$PWD/build -j 4

partitionByPhasedSNVs: PartitionByPhasedSNVs.cpp FastaIndex.h boost_1_66_0/stage/lib/libboost_program_options.a $(LIBBZ2)/libbz2.a $(LZMA)/liblzma.a  $(ZLIB)/build/lib/libz.a vcflib/lib/libvcflib.a
	$(CPP) -g $(CPPOPTS) $^ \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -I $(VCFLIB)/include -I $(HTSLIB) \
     -L $(BOOSTLIB) -l boost_program_options \
     -L $(VCFLIB)/lib -lvcflib  -L $(HTSLIB) -l hts -lpthread -L $(ZLIB)/build/lib -lz -L$(LZMA) -llzma -L$(LIBBZ2) -lbz2 -lpthread \
     -o $@ 


readToSNVList: ReadToSNVList.cpp PartitionTools.h FastaIndex.h SamUtils.h GenotypedRead.h SNVDB.h 
	$(CPP) $(CPPOPTS)  $< \
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
