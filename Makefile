all: vcflib/lib/libvcflib.a partitionByPhasedSNVs readToSNVList  

SEQAN=seqan/include
CONDA_LIB=testlib #$(CONDA_PREFIX)/lib
HTSINC=$(CONDA_PREFIX)/include
VCFLIB=vcflib

CPPOPTS= -O2
# -D_GLIBCXX_USE_CXX11_ABI=0


CPP=g++

#VCFLIB_INCLUDES := "-I $(abspath $(VCFLIB)/tabixpp/htslib)  -I$(VCFLIB)/include -I$(abspath vcflib/include) -L. -L$(CONDA_PREFIX)/lib -I$(abspath $(LIBBZ2)) -L$(abspath $(LIBBZ2))"

vcflib/lib/libvcflib.a:
	cd vcflib && make -j 8 libvcflib.a
	touch $@


#     -L $(VCFLIB)/lib -l vcflib \
#     -L $(CONDA_LIB) \
#     -l z -l pthread -lcurl -lssl -lcrypto -ldl  -llzma -lbz2 \

partitionByPhasedSNVs: PartitionByPhasedSNVs.cpp FastaIndex.h vcflib/lib/libvcflib.a
	$(CPP) $(CPPOPTS) -std=c++11 $<  \
     -o $@ \
     -I args \
     -I blasr \
     -I $(VCFLIB)/include -I $(HTSINC) \
     -L $(VCFLIB)/lib -l vcflib \
     -L $(VCFLIB)/tabixpp/htslib -l hts \
     -lpthread -lz -lm -llzma -lbz2 \
	   -lvcflib 

#partitionByPhasedSNVs: PartitionByPhasedSNVs.cpp FastaIndex.h vcflib/lib/libvcflib.a
#	$(CPP) $(CPPOPTS) -std=c++0x $<  \
#     -o $@ \
#     -I $(VCFLIB)/include -I $(HTSINC) \
#	   -L $(VCFLIB)/lib -l vcflib \
#     -L $(VCFLIB)/tabixpp/htslib -l hts \
#     -lpthread -lz -lm -llzma -lbz2 \
#	   -l vcflib 
#


readToSNVList: ReadToSNVList.cpp PartitionTools.h FastaIndex.h SamUtils.h GenotypedRead.h SNVDB.h 
	$(CPP) $(CPPOPTS)  $< \
     -I $(SEQAN) \
     -I args \
     -o $@ 

.PHONY: clean
clean:
	make -C $(LIBBZ2) -f Makefile-libbz2_so clean
	rm $(LIBBZ2)/libbz2.so
	make -C vcflib clean
	rm -rf partitionByPhasedSNVs readToSNVList
