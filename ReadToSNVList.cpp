#include "boost/program_options.hpp"
#include "seqan/vcf_io.h"
#include "seqan/seq_io.h"
#include "seqan/bam_io.h"
#include "seqan/align.h"
#include "seqan/basic.h"
#include <iostream>
#include <string>
#include <set>
#include "FASTASequence.h"
#include "FASTAReader.h"
#include "GenotypedRead.h"
#include "SamUtils.h"
#include "SampleTools.h"
#include "SNVDB.h"
namespace po = boost::program_options;

using namespace std;
using namespace seqan;
long min(long a, long b) {
	if (a <= b) {
		return a;
	} else {
		return b;
	}
}


typedef std::pair<int, GenotypedRead*> PosRead;


	

class CommandLineParser {
public:
	po::options_description desc;
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
	
	CommandLineParser() {
		vcfFileName = samFileName = refFileName = outFileName = nfqFileName = phasedSample = "";
		nfqOutFileName = "";
		minCoverage = 10;
		minFraction =0.25;
		refOffset=0;
		desc.add_options()
			("help", "Write help")
			("vcf", po::value<string>(), "VCF file. For now just het SNVs.")
			("nft", po::value<string>(), "Nucleotide frequency table.")
			("sam", po::value<string>(), "SAM file.")
			("ref", po::value<string>(), "Reference file.")
			("out", po::value<string>(), "Output file.")
			("nftOut", po::value<string>(), "Output filtered nucleotide frequency.")
			("minFraction", po::value<float>(), "Minimum fraction of coverge to represent.")
			("minCoverage", po::value<int>(), "Minimum absolute coverage for het call.")
			("phasedSample", po::value<string>(), "Sample to look up for phasing status.")
			("refOffset", po::value<int>(), "Offset into reference contig")
			;
	}

	template<typename T_Value>
	bool SetRequiredValue(string key, po::variables_map &vm, T_Value &value) {
		string res;
		if (vm.count(key.c_str())) {
			value =vm[key.c_str()].as<T_Value>(); 
			return true;
		} else {
			cout << desc << endl;
			cout << key << " is required.\n";
			exit(1);
		}
	}
	
	template<typename T_Value>
	void SetOptionalValue(string key, po::variables_map &vm, T_Value &value) {
		if (vm.count(key.c_str())) {
			value = vm[key.c_str()].as<T_Value>();
		}
	}		
		

	void ParseCommandLine(int ac, char* av[]) {
		po::variables_map vm;        
		try {
			po::store(po::parse_command_line(ac, av, desc), vm);
			po::notify(vm);    
			
			if (vm.count("help")) {
				cout << desc << "\n";
				return;
			}

			SetRequiredValue("sam", vm, samFileName);
			SetRequiredValue("ref", vm, refFileName);
			SetRequiredValue("out", vm, outFileName);
			SetOptionalValue("phasedSample", vm, phasedSample);
			SetOptionalValue("vcf", vm, vcfFileName);
			SetOptionalValue("nft", vm, nfqFileName);
			SetOptionalValue("nftOut", vm, nfqOutFileName);
			SetOptionalValue("minFraction", vm, minFraction);
			SetOptionalValue("minCoverage", vm, minCoverage);
			SetOptionalValue("refOffset", vm, refOffset);
		}

		catch(exception& e) {
			cerr << "error: " << e.what() << "\n";
			exit(1);
			return;
		}
		catch(...) {
			cerr << "Exception of unknown type!\n";
			exit(1);
		}
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

	po::options_description desc("Read to SNV list");

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
			vector<int>  snvPreBlock, snvPostBlock;
			vector<bool> retain;

			string alnChrom = genome[chromIndex].title;
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
					if (tAlnStr[i] != '-') {
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


							snvIndex = i;
							int allele = 2; // 0 = ref, 1 = alt, 2 = unknown

							if (chromSnps[curVar].nuc == qAlnStr[i]) {
								allele = 1;
							}
							else if (chromSnps[curVar].ref == qAlnStr[i]) {
								allele = 0;
							}
							else {
								allele = 2;
							}
							
							if (allele == 1) {
								chromSnps[curVar].nalt++;
							}
							else {
								chromSnps[curVar].nref++;
							}
							readSNVAlleles.push_back(allele);
						}
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
