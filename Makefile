all: partitionByPhasedSNVs readToSNVList 

CPPOPTS=  -g
CPP=g++ -std=c++14
ZLIB=zlib

partitionByPhasedSNVs: PartitionByPhasedSNVs.cpp FastaIndex.h 
	$(CPP) -g $(CPPOPTS) $^ \
     -lvcflib  -l hts -lz -llzma -lbz2 -lpthread \
     -o $@ 


readToSNVList: ReadToSNVList.cpp PartitionTools.h FastaIndex.h SamUtils.h GenotypedRead.h SNVDB.h 
	$(CPP) $(CPPOPTS)  $< \
     -I $(SEQAN) \
     -o $@ 

.PHONY: clean
clean:
	make -C $(LIBBZ2) -f Makefile-libbz2_so clean
	rm $(LIBBZ2)/libbz2.so
	make -C vcflib clean
	rm -rf partitionByPhasedSNVs readToSNVList
