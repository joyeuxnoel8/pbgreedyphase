#include "args.hxx"
#include <iostream>
#include <string>
#include <set>
#include "blasr/FASTASequence.h"
#include "blasr/FASTAReader.h"
#include "GenotypedRead.h"
#include "SamUtils.h"
#include "SampleTools.h"
#include "SNVDB.h"
#include "PartitionTools.h"

using namespace std;

long min(long a, long b) {
	if (a <= b) {
		return a;
	} else {
		return b;
	}
}



class CommandLineParser {
public:

	string vcfFileName;
	string samFileName;
	string refFileName;
	string outFileName;
	string nfqFileName;
	string nfqOutFileName;
	string phasedSample;
	float minFraction;
	int   minCoverage;
	int   refOffset;
	int minScoreDifference;
	bool  noRealign;
	CommandLineParser() {
		vcfFileName = samFileName = refFileName = outFileName = nfqFileName = phasedSample = "";
		nfqOutFileName = "";
		minCoverage = 10;
		minFraction =0.25;
		refOffset=0;
		minScoreDifference= 2;
		noRealign=false;

	}

	int ParseCommandLine(int ac, char* av[]) {
		args::ArgumentParser parser("Create a list of the SNVs a read overlaps", "");
    args::HelpFlag helpOpt(parser, "help", "Display this help menu", {'h', "help"});

    args::ValueFlag<string> samOpt(parser, "sam", "SAM file. ", {"sam"}, "", args::Options::Required);
    args::ValueFlag<string> refOpt(parser, "ref", "Reference. ", {"ref"}, "", args::Options::Required);
    args::ValueFlag<string> outOpt(parser, "out", "Output file. ", {"out"}, "", args::Options::Required);

    args::ValueFlag<string> vcfOpt(parser, "vcf", "VCF file. For now just het SNVs.", {"vcf"}, "");
    args::ValueFlag<string> nftOpt(parser, "nft", "Nucleotide frequency table. ", {"nft"}, "");
    args::ValueFlag<string> phasedSampleOpt(parser, "phasedSample", "Sample in vcf. ", {"sample"}, "");
		args::ValueFlag<int>	  minScoreDifferenceOpt(parser, 
																									"minScoreDifference",
																									"Minimum score difference between ref/alt realignment",
																									{"minScoreDifference"}, 2);
    args::Flag              noRealignOpt(parser, "no-realign", "Realign to target replaced by vcf difference.", {"no-realign"},false);
    
    args::ValueFlag<string> nftOutOpt(parser, "nftOut", "Nucleotide frequency table. ", {"nftOut"}, "");
		args::ValueFlag<float> minFractionOpt(parser, "minFraction", "Minimum fraction of coverge to represent.", {"minFraction"}, 0.25);
    args::ValueFlag<int> minCoverageOpt(parser, "minCoverage", "Minimum absolute coverage for het call.", {"minCoverage"},10);
    args::ValueFlag<int> refOffsetOpt(parser, "refOffset", "Offset into reference contig.", {"refOffset"}, 0);
		try {
			const std::vector<std::string> arguments(av + 1, av + ac);

			parser.ParseCLI(arguments);
    }
    catch (args::Completion e)
			{
        std::cout << e.what();
        return 0;
			}
    catch (args::Help)
			{
        std::cout << parser;
        return 0;
			}
    catch (args::ParseError e)
			{
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
			}

		vcfFileName = vcfOpt.Get();
		samFileName = samOpt.Get();
		refFileName = refOpt.Get();
		outFileName = outOpt.Get();
		nfqFileName = nftOpt.Get();
		nfqOutFileName = nftOutOpt.Get();
		phasedSample = phasedSampleOpt.Get();
		minFraction = minFractionOpt.Get();
		minCoverage = minCoverageOpt.Get();
		refOffset = refOffsetOpt.Get();
		minScoreDifference = minScoreDifferenceOpt.Get();
		noRealign = noRealignOpt.Get();
	 
	}
	
};

bool IsHetGenotype(const char* gt, int len){ 
	if (len < 3) {
		return false;
	}
	return (strncmp(gt, "0/1", 3) == 0 or 
					strncmp(gt, "1/0", 3) == 0 or
					strncmp(gt, "0|1", 3) == 0 or
					strncmp(gt, "1|0", 3) == 0);
}	

int StoreGenotype(const char* gt, int len, int &gt1, int &gt2) {
	// 
	// For now don't even try to handle the error.
	//
	if (len >= 3) {
		gt1 = (int) gt[0] - '0';
		gt2 = (int) gt[2] - '0';
		return 1;
	}
	else {
		return 0;
	}
}

bool OpenFile(string filename, ofstream &file) {
	if (filename == "") return false; 
	file.open(filename.c_str());
	if (!file.good()) {
		cout << "Could not open " << filename << endl;
		exit(1);
	}
	return true;
}


int MakeSNVDB(VcfFileIn &vcfFile, int sampleIndex, SNVDB &snvDb) {
	int nEntry = 0;
	VcfIOContext<>::TNameStore vcfContigNames;
	vcfContigNames = contigNames(context(vcfFile));

	while (!atEnd(vcfFile)) {
		VcfRecord vcfRecord;		
		VcfRecord record;
		readRecord(record, vcfFile);

		if (length(record.alt) == length(record.ref) && length(record.alt) == 1) {
			int recordLength = length(record.genotypeInfos[sampleIndex]);
			char *genotypeInfos = toCString(record.genotypeInfos[sampleIndex]);

			int genotypeInfosLength = length(record.genotypeInfos[sampleIndex]);
			if ( recordLength >= 3 and 
					 IsHetGenotype(genotypeInfos, genotypeInfosLength)) {
				int h1, h2;
				StoreGenotype(genotypeInfos, genotypeInfosLength, h1, h2);
				if (record.rID < length(vcfContigNames)) {
					// Add the alt snp to the database
					snvDb.AddSNV(toCString(vcfContigNames[record.rID]), record.beginPos, record.alt[0], record.ref[0], h1, h2);
				}
			}
			++nEntry;
		}
	}
	snvDb.Finalize();
	return nEntry;
}

int nucIndex[255];

int indexToNuc[4] = {'A','C','G','T'};


int main (int ac, char* av[]) {
	nucIndex[(int)'A'] = 0;
	nucIndex[(int)'C'] = 1;
	nucIndex[(int)'G'] = 2;
	nucIndex[(int)'T'] = 3;

	CommandLineParser	args;
	args.ParseCommandLine(ac, av);
	
	int sampleIndex;

	ifstream bamFileIn(args.samFileName.c_str());
	
	vector<FASTASequence> genome;
	FASTAReader reader;
	reader.SetToUpper();
	reader.Initialize(args.refFileName);
	reader.ReadAllSequences(genome);
		
	string regionChrom;
	int regionStart=0, regionEnd=-1;
	
	ofstream outFile(args.outFileName.c_str());

	//
	// Do some parsing of the vcf file.
	//
	SNVDB snvDb;
	if (args.vcfFileName != "") {
		VcfFileIn vcfFile(args.vcfFileName.c_str());
		VcfHeader vcfHeader;

		readHeader(vcfHeader, vcfFile);

		VcfIOContext<>::TNameStore vcfSampleNames;

		vcfSampleNames = sampleNames(context(vcfFile));
		int res = 0;
		if (args.phasedSample != "") {
			int res = LookupPhasedSampleIndex(vcfSampleNames, args.phasedSample, sampleIndex);	
			if (res == 0) {
				cout << "Did not find sample " << args.phasedSample << endl;
				exit(1);
			}
		}
		else {
			sampleIndex = 0;
		}


		int nEntry =  0;
		MakeSNVDB(vcfFile, sampleIndex, snvDb);
	}	


	//	readHeader(header, bamFileIn);
	SAMHeader samHeader;
	ReadHeader( bamFileIn, samHeader);

	map<string, int> refIndex;
	
	if (args.refOffset > 0) {
		FixedNameToIndexMap(samHeader, genome, refIndex);
	}
	else {
		BuildNameToIndexMap(samHeader, genome, refIndex);
	}
	if (refIndex.size() == 0) {
		cout << "ERROR, there is no reference map. Perhaps \"samtools view\" was used without -h?" << endl;
		exit(0);
	}

	
	if (args.nfqFileName != "") {
		ifstream nfq(args.nfqFileName.c_str());
		ofstream nfqOut;
		if (args.nfqOutFileName != "") {
			nfqOut.open(args.nfqOutFileName.c_str());
		}
		string contig;
		int pos;
		int nuc[6];
		int ndel,nins;
		
		string insStr, ambigStr;
		int nfqIndex = 0;
		string line;
		while (getline(nfq, line)) {
			stringstream nfqStrm(line);
			if (!(nfqStrm >> contig >> pos >> nuc[0] >> nuc[1] >> nuc[2] >> nuc[3] >> nuc[4] >> nuc[5]) ) {
				break;
			}
			pos -= args.refOffset + 1;
			//
			//  Find the SNV.
			//
			int total = 0;
			int hasMin[4] = {0,0,0,0};
			int hasFrac[4] = {0,0,0,0};
			int i;
			total = 0;
			for (i = 0; i < 6;i++) {
				total+= nuc[i];
			}
			float f = total*args.minFraction;
			int nFrac = 0;
			int nMin = 0;
			int maxCoverage = 0;
			int maxCoverageIndex = -1;
			for (i =0; i < 4;i++) {
				if (nuc[i] >= f) {
					nFrac+=1;
					hasFrac[i] = 1;
				}
				if (nuc[i] >= args.minCoverage) {
					nMin +=1;
					hasMin[i] = 1;
				}
				if (nuc[i] > maxCoverage) {
					maxCoverage = nuc[i];
					maxCoverageIndex = i;
				}
			}
			++nfqIndex;
			//			if ((args.minFraction > 0 and nFrac != 2) or 
			//					nMin != 2) {
			//			cerr << nMin<< endl;
			if (nMin != 2) {
				cerr << "Skipping line: " << nfqIndex << " " << nFrac << " " << nMin << " " << contig << " " << pos << " " << nuc[0] << " " << nuc[1] << " " << nuc[2] <<" " << nuc[3] << " " << nuc[4] << " " << nuc[5] << endl;
				continue;
			}
			if (args.nfqOutFileName != "") {
				nfqOut << line << endl;
			}
			char refNuc;
			int chromIndex = refIndex[contig];
			refNuc = genome[chromIndex].seq[pos];
			int refIndex = nucIndex[refNuc];
			char queryNuc = '\0';
			for (i = 0; i < 4; i++) {
				if (nuc[i] >= args.minCoverage and
						i != refIndex) {
					queryNuc = indexToNuc[i];
				}
			}
			assert(queryNuc != '\0');
			
			snvDb.AddSNV(contig, pos, queryNuc, refNuc);
		}
		cerr << "ended on  " << pos << endl;
		snvDb.Finalize();
	}
	cerr << "SNVDB " << snvDb.size << endl;

	//
	// Do some parsing of the sam file.
	//

	// Copy header.
	typedef Row<Align<IupacString> >::Type TRow;                 // gapped sequence type

	
	int maybeStored = 0;
	while (bamFileIn)  {
		SAMRecord samRecord;
		bool result;
		result = ReadRecord(bamFileIn, samRecord);
		if (result == false) {
			break;
		}
		Align<IupacString> align;
		
		if (samRecord.chrom != "*") {
			samRecord.refPos = samRecord.refPos - args.refOffset;

			int firstOp = 0;
			int totalClipped = 0;
			int refAlnLength = samRecord.GetRefAlignLength();
			int refAlnStart  = samRecord.refPos;
			int chromIndex = refIndex[samRecord.chrom];
			if (chromIndex < 0) {
				continue;
			}
			char *refChrom = (char*) genome[chromIndex].seq;
			string tAlnStr, qAlnStr;
			MakeAlignStrings(samRecord, refChrom, regionStart, regionEnd, qAlnStr, tAlnStr);
			int l1 = tAlnStr.size();
			int l2 = qAlnStr.size();
			
			int minLength = min(l1, l2);
			
			int qPos = 0;
			int tPos = 0;
			int i;
			int nSNP = 0;

			vector<int>  chromSNVPos;
			vector<char> chromSNVChar;
			vector<char> chromSNVRef;
			vector<char> readSNVChar;
			vector<int>  readSNVAlleles;
			vector<int>  snvPreBlock, snvPostBlock,alnScoreDiff;
			vector<bool> retain;

			string alnChrom = genome[chromIndex].GetName();
			int startVar, endVar, curVar;
			bool foundBounds;
			//
			// Save start var and end var, the first and last variants that overlap
			// the alignment of this read.
			//
			
			foundBounds = snvDb.QueryBounds(alnChrom, refAlnStart, refAlnStart + refAlnLength, startVar, endVar);
			bool setPostBlock = false;
			int snvIndex = 0;
			GenotypedRead *read = new GenotypedRead;
			read->chrom    = alnChrom;
			read->samLine  = samRecord.samLine;
			read->alnStart = refAlnStart;
			read->alnEnd   = refAlnStart + refAlnLength;
			int sitesGenotyped = 0;
			int sitesUnknown    = 0;
			int span = endVar - startVar;
			
			if (foundBounds and startVar != endVar) {
				curVar = startVar;
				SNVDB::SNVs &chromSnps = snvDb.db[alnChrom];
				int lastGap = 0;
				//
				// Scan the alignment for positions that overlap SNVs.
				//
				for (i = 0; i < minLength && curVar < endVar; i++) {
					//
					// Bookkeeping to tabulate accuracy of this read.
					//
					if (tAlnStr[i] != qAlnStr[i]) {
						++nSNP;
					}

					//
					// Make sure checking the most up to date position of the snv.
					//
					while (curVar < endVar and chromSnps[curVar].pos < tPos + refAlnStart) {
						curVar ++;
					}
						
					//
					// This is the condition that a SNV was found.
					//
					if (curVar < endVar and chromSnps[curVar].pos == tPos + refAlnStart) {
						chromSNVPos.push_back(tPos + refAlnStart);
						chromSNVRef.push_back(genome[chromIndex].seq[tPos + refAlnStart]);
						chromSNVChar.push_back(chromSnps[curVar].nuc);
						readSNVChar.push_back(qAlnStr[i]);

						//
						// Store how long of match was before the mismatch, and
						// signal to store how long after.
						//

						int p = i-1;
						while (p > 0 and tAlnStr[p] != '-' and qAlnStr[p] != '-') { p--; }
						snvPreBlock.push_back(i-p-1);
						p = i+1;
						while (p < tAlnStr.size() and tAlnStr[p] != '-' and qAlnStr[p] != '-') { p++; }
						snvPostBlock.push_back(p-i-1);
						int w = 10;
						string tPre, tSuf, qPre, qSuf;
						int refScore =0, altScore=0;
						snvIndex = i;
						int allele = 2; // 0 = ref, 1 = alt, 2 = unknown
							
						if (args.noRealign == false) {
							if (GetUngappedPrefix(tAlnStr, i, w, tPre) and
									GetUngappedPrefix(qAlnStr, i, w, qPre) and
									GetUngappedSuffix(tAlnStr, i, w, tSuf) and
									GetUngappedSuffix(qAlnStr, i, w, qSuf)) {
								string tRef = tPre + chromSnps[curVar].ref + tSuf;
								string tAlt = tPre + chromSnps[curVar].nuc + tSuf;
								string qStr = qPre + qAlnStr[i] + qSuf;
									
									
								string refTStr, refQStr, altTStr, altQStr;
								refScore=SWAlign(tRef, qStr, -2,1,2, refQStr, refTStr);
								altScore=SWAlign(tAlt, qStr, -2,1,2, altQStr, altTStr);
									
								/*
								//
								// Lots of debugging information
								//
								cout << "REF\t \t\tALT" << endl;
								cout << tPre << ' ' << chromSnps[curVar].ref << ' ' << tSuf << "\t"
								<< tPre << ' ' << chromSnps[curVar].nuc << ' ' << tSuf << endl;
								cout << qPre << ' ' << qAlnStr[i] << ' ' << qSuf << "\t"
								<< qPre << ' ' << qAlnStr[i] << ' ' << qSuf<< endl;
								cout << ungappedPrefix << "," << ungappedSuffix << endl;
								cout <<  refScore << endl;
								cout << refQStr << endl;
								cout << refTStr << endl << endl;
								cout << altScore << endl;
								cout << altQStr << endl;
								cout << altTStr << endl << endl;
								*/
								
							}

							int scoreDiff = refScore-altScore;
							alnScoreDiff.push_back(scoreDiff);

							//							if (chromSnps[curVar].nuc == qAlnStr[i]) {
							if (scoreDiff < -args.minScoreDifference) {
								allele = 1;
							}
							//							else if (chromSnps[curVar].ref == qAlnStr[i]) {
							else if (scoreDiff > args.minScoreDifference) {
								allele = 0;
							}
							else {
								allele = 2;
							}
						}
						else {
							//
							// assign allele directly from the pileup
							//
								
							if (chromSnps[curVar].nuc == qAlnStr[i]) {
								allele = 1;
							}
							else if (chromSnps[curVar].ref == qAlnStr[i]) {
								allele = 0;
							}
							else {
								allele = 2;
							}
							refScore = altScore;
						}
							
						if (allele == 1) {
							chromSnps[curVar].nalt++;
						}
						else {
							chromSnps[curVar].nref++;
						}
						readSNVAlleles.push_back(allele);
					}

					if (tAlnStr[i] == '-' or qAlnStr[i] == '-') {
						lastGap = i;
					}
					if (tAlnStr[i] != '-') {
						tPos ++;
					}
					if (qAlnStr[i] != '-') {
						qPos ++;
					}
				}
				
				outFile << samRecord.title << "\t" << chromSNVPos.size() << "\t" << samRecord.chrom;
				for (i = 0; i < chromSNVPos.size(); i++) {
					outFile << "\t";
					outFile << chromSNVPos[i] << "," 
									<< chromSNVRef[i] << ","
									<< chromSNVChar[i] << ","
									<< readSNVChar[i] << ","
									<< snvPreBlock[i] << ","
									<< snvPostBlock[i];
				}
				outFile << endl;

			}
		}
	}

	return 0;
}
