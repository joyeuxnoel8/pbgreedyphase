#include "boost/program_options.hpp"
#include "seqan/vcf_io.h"
#include "seqan/seq_io.h"
#include "seqan/bam_io.h"
#include "seqan/align.h"
#include "seqan/basic.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <set>
#include "FASTASequence.h"
#include "FASTAReader.h"
#include "SampleTools.h"
#include "FastaIndex.h"
#include "PartitionTools.h"

using namespace std;

#include "Variant.h"
using namespace seqan;
namespace po = boost::program_options;


long min(long a, long b) {
	if (a <= b) {
		return a;
	} else {
		return b;
	}
}

class SAMHeader {
public:
	vector<string> lines;
	vector<string> chroms;
	vector<int> lengths;
};

class GenotypedRead {
public:
	string chrom;
	int alnStart, alnEnd;
	vector<int> pos;
	vector<char> genotype; 
	//
	// Helpful for weighing snv by accuracy
	//
	vector<int> preBlock;
	vector<int> postBlock;
	string samLine;
	void PrintHaplotype() {
		int i;
		for (i =0;i<genotype.size();i++) {
			cout << (int) genotype[i] <<",";
		}
		cout << endl;
	}
};

typedef std::pair<int, GenotypedRead*> PosRead;

void ReadHeader(istream &in, SAMHeader &header) {
	while (in.peek() == '@') {
		string line;
		getline(in, line);
		if (line == "") {
			return;
		}
		header.lines.push_back(line);
		stringstream tagStrm(line);
		string tagKey, tagVal;
		string headerType;
		tagStrm >> headerType;
		if (headerType != "@SQ") {
			continue;
		}
		tagStrm.get();
		tagKey.resize(2);
		string sn("SN");
		string ln("LN");

		while (tagStrm) {
			string tagKV;
			if ( !(tagStrm >> tagKV) ) break;
			stringstream tagKVStrm(tagKV);
			tagKVStrm.get(&tagKey[0], 3);
			tagKVStrm.get();
			tagKVStrm >> tagVal;
			
			if (tagKey == sn) {
				header.chroms.push_back(tagVal);
			}
			if (tagKey == ln) {
				header.lengths.push_back(atoi(tagVal.c_str()));
			}
		}
	}
}

void PrintSAMHeader(SAMHeader &header, ostream &out) {
	int i;
	for (i = 0; i < header.lines.size(); i++) {
		out << header.lines[i] << endl;
	}
}


void BuildNameToIndexMap(SAMHeader &header, vector<FASTASequence> &genome, map<string, int> &refIndex) {
	int i;
	int j;

	for (i = 0; i < header.chroms.size(); i++) {
		for (j = 0; j < genome.size(); j++) {
			if (header.chroms[i] == genome[j].title) {
				refIndex[header.chroms[i]] = j;
				break;
			}
		}
		if (j == genome.size()) {
			refIndex[header.chroms[i]] = -1;
		}
	}
}

bool IsMatch(char c) {
	return c == 'M' or c == '=' or c == 'X';
}

class SAMRecord {
public:
	string title;
	int flag;
	int refPos;
	int mapqv;
	string cigar, seq;
	string chrom;
	vector<int> lengths;
	vector<char> ops;
	string samLine;
	int tLen;
	int GetRefAlignLength() {
		int refAlignLength=0, i;
		for (i = 0; i < lengths.size(); i++) {
			if (IsMatch(ops[i]) or ops[i] == 'D') {
				refAlignLength += lengths[i];
			}
		}
		assert(tLen == refAlignLength);
		return refAlignLength;
	}
};
	
bool ReadRecord(istream &in, SAMRecord &record) {
	int intdummy;
	string dummy;
	string line;
	getline(in, line);
	if (line == "") 
		return false;
	stringstream lineStrm(line);
	record.samLine = line;

	lineStrm >> record.title
					 >> record.flag 
					 >> record.chrom 
					 >> record.refPos 
					 >> record.mapqv 
					 >> record.cigar 
					 >> dummy 
					 >> intdummy 
					 >> record.tLen
					 >> record.seq;
	// sam is 1 based
	record.refPos -=1;
	stringstream cigarstrm(record.cigar);
	while (cigarstrm) {
		int length;
		char op;
		if ( !( cigarstrm >> length) ) break;
		cigarstrm.get(op);
		if (op == EOF) {
				break;
		}
		record.lengths.push_back(length);
		record.ops.push_back(op);
	}
	return true;
}

int MakeAlignStrings(SAMRecord &record, char *ref, int regionStart, int regionEnd, string &queryAln, string &refAln) {
	//
	// First compute the length of the output strings.
	//
	int i;
	int alnLength = 0;
	for (i = 0; i < record.lengths.size(); i++) {
		if (record.ops[i] != 'H' and record.ops[i] != 'S') {
			alnLength += record.lengths[i];
		}
	}
	queryAln.resize(alnLength);
	refAln.resize(alnLength);
	//
	// Now store output strings.
	//
	char *query = &record.seq[0];
	int queryPos = 0;
	int refPos = record.refPos;
	int queryAlnPos = 0;
	int refAlnPos = 0;
	for (i = 0; i < record.lengths.size(); i++) {
		if (record.ops[i] == 'H') {
			continue;
		}
		if (record.ops[i] == 'S') {
			queryPos += record.lengths[i];
		}
		else {
			break;
		}
	}
	//
	// Advance pointers through the reference 
	for (; refPos < regionStart && i < record.lengths.size(); refPos < regionStart && i++) {
		if (record.ops[i] == 'M' or record.ops[i] == 'X' or record.ops[i] == '=') {
			int copyLength = min(regionStart - refPos, record.lengths[i]);
			refPos   += copyLength;
			queryPos += copyLength;
			record.lengths[i] -= copyLength;
		}
		if (record.ops[i] == 'D') {
			int copyLength = min(regionStart - refPos, record.lengths[i]);			
			refPos += copyLength;
			record.lengths[i] -= copyLength;
		}
		if (record.ops[i] == 'I') {
			queryPos += record.lengths[i];
		}
	}
	int startRefPos = refPos;
	
	for (;i < record.lengths.size() && refPos < regionEnd ; i++) {
		if (record.ops[i] == 'H') {
			continue;
		}
		if (record.ops[i] == 'S') {
			continue;
		}
		if (record.ops[i] == 'M' or record.ops[i] == 'X' or record.ops[i] == '=') {
			if (regionEnd != -1 and refPos - startRefPos > regionEnd - regionStart) {
				queryAln = "";
				refAln = "";
				return 0;
			}
			int copyLength = min(record.lengths[i], regionEnd - refPos);
			strncpy(&queryAln[queryAlnPos], &query[queryPos], copyLength);
			strncpy(&refAln[refAlnPos], &ref[refPos-regionStart], copyLength);
			queryAlnPos += copyLength;
			queryPos    += copyLength;
			refAlnPos   += copyLength;
			refPos      += copyLength;
		}
		if (record.ops[i] == 'I') {

			strncpy(&queryAln[queryAlnPos], &query[queryPos], record.lengths[i]);
			queryAlnPos += record.lengths[i];
			queryPos += record.lengths[i];
			int j;
			for (j = 0; j < record.lengths[i]; j++, refAlnPos++) {
				refAln[refAlnPos] = '-';
			}
		}
		if (record.ops[i] == 'D') {
			int copyLength = min(record.lengths[i], regionEnd - refPos);
			/*
			if (regionEnd != -1 and refPos - record.refPos > regionEnd - regionStart) {
				queryAln = "";
				refAln = "";
				return 0;
				}*/
			//			cout << i << " " << record.lengths[i] << " " << record.lengths.size() << endl;
			assert(refAlnPos >= 0);
			assert(refPos-regionStart >= 0);
			
			strncpy(&refAln[refAlnPos], &ref[refPos-regionStart], copyLength);
			refAlnPos += copyLength;
			refPos += copyLength;
			int j;
			for (j = 0; j < copyLength; j++, queryAlnPos++) {
				queryAln[queryAlnPos] = '-';
			}
		}
	}
	for (int i = 0; i < refAln.size(); i++) {
		refAln[i] = toupper(refAln[i]);
	}
	return 1;
}


bool ParseRegion(string &region, string &chrom, int &start,int &end){
	string noCommas = "";
	int i;
	for (i = 0; i < region.size(); i++) {
		if (region[i] != ',') {
			noCommas.push_back(region[i]);
		}
	}
	

	int cpos =noCommas.find(":");
	if (cpos == noCommas.npos) {
		return false;
	}
	else {
		chrom = noCommas.substr(0,cpos);
		cerr << "rgn: " << chrom << endl;
	}
	cpos++;
	stringstream posStrm(noCommas.substr(cpos));
	if ((posStrm >> start).eof()) {
		return false;
	}
	cerr << "rgn start: " << start << endl;	
	posStrm.get();
	if ((posStrm >> end).bad()) {
		return false;
	}
	return true;
}
	

class CommandLineParser {
public:
	po::options_description desc;
	string vcfFileName;
	string samFileName;
	string refFileName;
	string h1FileName;
	string h2FileName;
	string sample;
	string phaseTag;
	int verbosity;
	int block;
	int minGenotyped;
	int maxUnknown;
	int minDifference;
	int minScoreDifference;
	string unassigned;
	string region;
	string summaryFile;
	string phaseStatsFileName;
	string chromosome;
	bool assumeAutozygous;
	int padding;
	int nwWindow;
	int writeInterval;
	CommandLineParser() {
		h1FileName = "";
		h2FileName = "";
		verbosity = 0;
		block = 3;
		minGenotyped = 2;
		minScoreDifference= 1;
		maxUnknown = 6;
		unassigned = "";
		region = "";
		summaryFile = "";
		phaseStatsFileName = "";
		minDifference = 0;
		nwWindow = 10;
		phaseTag = "";
		padding=0;
		vcfFileName="";
		refFileName="";
		samFileName="";
		int    regionPadding;
		writeInterval=0;
		assumeAutozygous=false;
		desc.add_options()
			("help", "Write help")
			("vcf", po::value<string>(), "VCF file. For now just het SNVs.")
			("sam", po::value<string>(), "SAM file")
			("ref", po::value<string>(), "Reference file")
			("tag", po::value<string>(), "Store haplotype in SAM optional field with this 2-character tag.")
			("h1", po::value<string>(), "Haplotype 1")
			("h2", po::value<string>(), "Haplotype 2")
			("verbosity", po::value<int>(),"Verbosity")
			("pad", po::value<int>(), "padding around region.")
			("block", po::value<int>(), "Minimum adjacent block size to record SNV status")
			("minDifference", po::value<int>(), "Minimum difference between genotype count")
			("minScoreDifference", po::value<int>(), "Minimum score difference between ref/alt realignment")
			("nw-window", po::value<int>(), "Window to do Needleman-Wunsch alignment in")			
			("sample", po::value<string>(), "Sample to look up for phasing status.")
			("maxUnknown", po::value<int>(), "Maximum sites with unknown/alt genotype")
			("minGenotyped", po::value<int>(), "Minimum genotyped sites per read")
			("unassigned", po::value<string>(), "Output unassigned reads here")
			("rgn", po::value<string>(), "Region of reference to phase")
			("summary",  po::value<string>(), "Write a summary of phased reaads to this file.")
			("phaseStats",  po::value<string>(), "Write the number of h0/h1 matches to this file per read.")
			("null-assumes-autozygous",  "If there are no values in the VCF, assume autozygous region.")
			("writeInterval", "Write an interval instead of the whole sam line.")
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

	void SetOptionalFlag(string key, po::variables_map &vm, bool &value) {
		if (vm.count(key.c_str())) {
   		value = true;
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
			if (vm.count("writeInterval")) {
				writeInterval = 1;
			}
			SetRequiredValue("vcf", vm, vcfFileName);
			SetRequiredValue("sam", vm, samFileName);
			SetRequiredValue("ref", vm, refFileName);

			
			SetOptionalValue("sample", vm, sample);
			//
			// Parse optional sargs.
			//
			SetOptionalValue("pad", vm, padding);			
			SetOptionalValue("rgn", vm, region);
			SetOptionalValue("tag", vm, phaseTag);
			SetOptionalValue("h1", vm, h1FileName);
			SetOptionalValue("h2", vm, h2FileName);
			SetOptionalValue("verbosity", vm, verbosity);
			SetOptionalValue("block", vm, block);
			SetOptionalValue("minDifference", vm, minDifference);			
			SetOptionalValue("maxUnknown", vm, maxUnknown);
			SetOptionalValue("minGenotyped", vm, minGenotyped);
			SetOptionalValue("unassigned", vm, unassigned);
			SetOptionalValue("phaseStats",vm, phaseStatsFileName);
			SetOptionalValue("chrom",vm, chromosome);
			SetOptionalValue("summary",vm, summaryFile);
			SetOptionalValue("nw-window",vm,nwWindow);
			SetOptionalFlag("null-assumes-autozygous",vm, assumeAutozygous);
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



class SNV {
public:
	int pos;
	char nuc, ref;
	int h1, h2;
	SNV(int p) {
		pos = p;
		nuc = '\0';
		ref = '\0';
		h1 = h2 = 0;
	}

	SNV(int p, char n, char r, int _h1, int _h2) {
		pos = p;
		nuc = n;
		ref = r;
		h1  = _h1;
		h2  = _h2;
	}

	int operator<(int rhs) const {
		return pos < rhs;
	}

	int operator<(const SNV &rhs) const {
		return pos < rhs.pos;
	}
	int operator==(const SNV &rhs) const {
		return pos == rhs.pos;
	}
	SNV& operator=(const SNV &rhs) {
		pos = rhs.pos;
		nuc = rhs.nuc;
		ref = rhs.ref;
		h1  = rhs.h1;
		h2  = rhs.h2;
		return *this;
	}
	
	bool GetPhase(int allele, int &phase) {
		if (allele == h1) { phase = 0; return true; }
		if (allele == h2) { phase = 1; return true; }
		phase = 2;
		return false;
	}
};

bool IsPhasedHetGenotype(string &gt) {
	if (gt.size() < 3) {
		return false;
	}
	return (gt == "0|1" or gt == "1|0");
}	

int StoreGenotype(string gt, int &gt1, int &gt2) {
	// 
	// For now don't even try to handle the error.
	//
	if (gt.size() == 3) {
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

class SNVDB {
public:
	
	int size;
	typedef vector<SNV> SNVs;
	typedef map<string, SNVs > ChromDB;
	ChromDB db;

	
	SNVDB() {
		size = 0;
	}

	string MakeKey(string chrom, int pos, char nuc) {
		stringstream snvstrm;
		snvstrm << chrom << "_" << pos << "_" << nuc;
		return snvstrm.str();
	}

	
	void AddSNV(string chrom, int pos, char nuc, char ref, int h1, int h2) {
		
		if (db.find(chrom) == db.end()) {
			db[chrom] = SNVs();
		}

		db[chrom].push_back(SNV(pos, nuc, ref, h1, h2));
		size++;
		if (size % 100000 == 0) {
			cerr << "snvdb: " << size << endl;
		}
	}
	
	bool QuerySNV(string chrom, int pos, char nuc) {
		cout << "NOT WRITTEN";
		assert(0);
	}

	void Finalize() {
		ChromDB::iterator dbIt;
		for (dbIt = db.begin(); dbIt != db.end(); ++dbIt) {

			std::sort( dbIt->second.begin(), dbIt->second.end() );
		}
	}

	bool QueryBounds(string chrom, int start, int end,  int &firstVar, int &endVar) {
		if (db.find(chrom) == db.end()) {
			firstVar = endVar = 0;
			return false;
		}
		// Get quick access to the vector
		SNVs &vect = db[chrom];
		// it doesn't matter what nuc this is.
		SNV startQuery(start);
		SNV endQuery(end);
		SNVs::iterator firstIt = std::lower_bound(vect.begin(), vect.end(), start);
		SNVs::iterator lastIt = std::lower_bound(vect.begin(), vect.end(), end);
		while (lastIt != vect.end() and (*lastIt).pos < end) {
			++lastIt;
		}
		firstVar = firstIt - vect.begin();
		endVar   = lastIt  - vect.begin();
		return true;
	}

};


bool AnyHaplotypeOverlaps(int overlaps[],int len) {
	int i;
	for (i = 0; i < len; i++) {
		if (overlaps[i]) { return true;}
	}
	return false;
}

template<typename T>
void FilterList(vector<T> &v, vector<bool> &retain) {
	int i, c;
	c = 0;
	for (i = 0; i < retain.size(); i++) {
		if (retain[i] == true) {
			v[c] = v[i];
			c++;
		}
	}
}

int main (int ac, char* av[]) {
	po::options_description desc("Partition by SNV options");

	CommandLineParser	args;
	args.ParseCommandLine(ac, av);
	int sampleIndex;

	ifstream bamFileIn(args.samFileName.c_str());
	
	vector<FASTASequence> genome;
	FASTAReader reader;
	reader.SetToUpper();
	
	GFFastaIndex fastaIndex;
	fastaIndex.Initialize(args.refFileName);
		
	string regionChrom;
	int regionStart=0, regionEnd=-1;
	if (args.region != "") {
		bool res;
		res = ParseRegion(args.region, regionChrom, regionStart, regionEnd);
		cerr << res << endl;
		if (!res) {
			cerr << "Could not parse " << args.region << endl;
			exit(1);
		}
		regionStart -= args.padding;
		regionEnd   += args.padding;
	}
	else {
		reader.Initialize(args.refFileName);
		reader.ReadAllSequences(genome);
	}

	
	sampleIndex = 0;
	ofstream haplotypeOut[2];
	ofstream phaseStatsOut;
	if (args.phaseStatsFileName != "") {
		phaseStatsOut.open(args.phaseStatsFileName.c_str());
	}
	int nEntry =  0;
	SNVDB snvDb;

	bool useH1 = OpenFile(args.h1FileName, haplotypeOut[0]);
	bool useH2 = OpenFile(args.h2FileName, haplotypeOut[1]);

	vcflib::VariantCallFile vcf;
	vcf.open(args.vcfFileName);
	vcflib::Variant record(vcf);
	//	while (!atEnd(vcfFile)) {
	int i;
	while (vcf.getNextVariant(record)) {

		string ref = record.ref;
		string alt = record.alt[0];
		string chrom = record.sequenceName;
		int pos = record.position;
			

		if (ref.size() == alt.size() and ref.size() == 1) {
			int recordLength = ref.size();
			map<string, vector<string> > *sample;
			if (args.sample == "") {
				sample = &record.samples.begin()->second;
			}
			else {
				//
				// Do some work to find the sample
				//
				vcflib::Samples::iterator sampIt = record.samples.find(args.sample);
				if (sampIt == record.samples.end()) {
					cout << "ERROR, sample " << args.sample << " was not in the vcf file." << endl;
					int s;
					cout << "Samples include: ";
					for (s = 0; s < record.sampleNames.size(); s++) {
						cout << record.sampleNames[s] << ",";
					}
					cout << endl;
					exit(1);
				}
				sample = &sampIt->second;
			}
			if (sample->find("GT") != sample->end()) {
				string gt = (*sample)["GT"][0];
				if (gt.size() == 3 and
						IsPhasedHetGenotype(gt) ) {
					int h1, h2;
					StoreGenotype(gt,h1, h2);
					snvDb.AddSNV(chrom, pos-1, ref[0], alt[0], h1, h2);
				}
			}

			++nEntry;
		}
	}
	snvDb.Finalize();
	
	cerr << "SNVDB " << snvDb.size << endl;

	//
	// Do some parsing of the sam file.
	//

	// Copy header.

	//	readHeader(header, bamFileIn);
	SAMHeader samHeader;
	ReadHeader( bamFileIn, samHeader);

	map<string, int> refIndex;
	BuildNameToIndexMap(samHeader, genome, refIndex);
	if (refIndex.size() == 0) {
		cout << "ERROR, there is no reference map. Perhaps \"samtools view\" was used without -h?" << endl;
		exit(0);
	}
	ofstream unassigned;
	if (args.writeInterval == false){ 
		for (i = 0; i< (useH1 + useH2); i++) {
			PrintSAMHeader(samHeader, haplotypeOut[i]);
		}
		unassigned.open(args.unassigned.c_str());
		if (args.unassigned != "") {
			PrintSAMHeader(samHeader, unassigned);
		}
	}

	typedef Row<Align<IupacString> >::Type TRow; // gapped sequence type

	int maybeStored = 0;

	string refChrom, prevRef;
	prevRef = "";

	int nH1=0,nH2=0,nUnk=0;
	int bH1=0,bH2=0,bUnk=0;
	bool writeField = false;
	bool splitHaplotypes = false;

	if (args.phaseTag != "") {
		writeField = true;
		splitHaplotypes = false;
	}
	if (args.h1FileName != "" and args.h2FileName != "") {
		splitHaplotypes = true;
	}

	int refStart = 0;
	int refRegionStart = 0;
	if (regionStart < 0) {
			regionStart = 0;
	}
	
	if (args.region != "" &&
			fastaIndex.fai.find(regionChrom) != fastaIndex.fai.end()) {
		int contigLength = fastaIndex.fai[regionChrom][0];
		if (contigLength < regionEnd || regionEnd == -1) {
			regionEnd = contigLength;
			cout << "set region end." <<endl;
			exit(0);
		}
		fastaIndex.GetSeq(refChrom, regionChrom, regionStart, regionEnd);
	}
	
	int alnIndex = 0;
	while (bamFileIn)  {
		SAMRecord samRecord;
		bool result;
		result = ReadRecord(bamFileIn, samRecord);
		alnIndex +=1;
		
		if (result == false) {

			break;
		}
		
		Align<IupacString> align;
		
		if (samRecord.chrom != "*") {
			int firstOp = 0;
			int totalClipped = 0;
			int refAlnLength = samRecord.GetRefAlignLength();
			int refAlnStart  = max(samRecord.refPos, regionStart);
			int chromIndex   = refIndex[samRecord.chrom];
			
			string tAlnStr, qAlnStr;
			char *refPtr;
			if (args.region == "") {
				refPtr = (char*) genome[chromIndex].seq;
				regionStart = 0;
				regionEnd=genome[chromIndex].length;
			}
			else {
				refPtr = &refChrom[0];
			}
				
			int res = MakeAlignStrings(samRecord, refPtr, regionStart, regionEnd, qAlnStr, tAlnStr );
			if (res == 0) {
				cerr << "Could not make aln strings " << refPtr << " " << regionStart << " " << regionEnd << endl;
				continue;
			}
			
			int l1 = tAlnStr.size();
			int l2 = qAlnStr.size();
			
			int minLength = min(l1, l2);
			
			int qPos = 0;
			int tPos = 0;
			int i;

			string alnChrom = samRecord.chrom;
			int startVar, endVar, curVar;
			bool foundBounds;
			foundBounds = snvDb.QueryBounds(alnChrom, refAlnStart, refAlnStart + refAlnLength, startVar, endVar);
			bool setPostBlock = false;
			int snvIndex = 0;
			/*
				GenotypedRead *read = new GenotypedRead;
				read->chrom = alnChrom;
				read->samLine = samRecord.samLine;
				read->alnStart = refAlnStart;
				read->alnEnd   = refAlnStart + refAlnLength;
			*/
			int sitesGenotyped = 0;
			int sitesUnknown    = 0;
			int span = endVar - startVar;

			int haplotype = 2;
			int h0=0,h1=0,h2=0;
						
			if (args.verbosity > 0) {
				cerr << samRecord.title << " overlaps " << endVar - startVar << " variants" << endl;
			}
			if (foundBounds and startVar != endVar) {
				
				curVar = startVar;
				SNVDB::SNVs &chromSnps = snvDb.db[alnChrom];
				
				int lastGap = 0;
				for (i = 0; i < minLength && curVar < endVar; i++) {

					while (curVar < endVar and chromSnps[curVar].pos < tPos + refAlnStart) {
						curVar ++;
					}
					  
					if (curVar < endVar and chromSnps[curVar].pos == tPos + refAlnStart) {

						string tPre, tSuf, qPre, qSuf;
						int refScore =0, altScore=0;
						string refTStr, refQStr, altTStr, altQStr;
						string tRef, tAlt, qStr;
						int t;
						if (GetUngappedPrefix(tAlnStr, i, args.nwWindow, tPre) and
								GetUngappedPrefix(qAlnStr, i, args.nwWindow, qPre) and
								GetUngappedSuffix(tAlnStr, i, args.nwWindow, tSuf) and
								GetUngappedSuffix(qAlnStr, i, args.nwWindow, qSuf)) {
							tRef = JoinPrefixSuffix(tPre, chromSnps[curVar].ref, tSuf);
							tAlt = JoinPrefixSuffix(tPre,chromSnps[curVar].nuc,tSuf);
							qStr = JoinPrefixSuffix(qPre,qAlnStr[i],qSuf);
							
							if (args.verbosity > 1) {
								cerr << tRef << endl
										 << tAlt << endl
										 << qStr << endl;
							}
							
							refScore=SWAlign(tRef, qStr, -2, 1, 4, refQStr, refTStr);
							altScore=SWAlign(tAlt, qStr, -2, 1, 4, altQStr, altTStr);
							if (args.verbosity > 1) {
								cerr << refQStr << " - " << altQStr << endl
										 << refTStr << " - " << altTStr << endl;
							}
							
						} else {
							if (args.verbosity > 1) {
								cerr << "Could not run overlap." << endl;
							}
						}
							
						int scoreDiff = refScore-altScore;

							int allele = 2; // 0 = ref, 1 = alt, 2 = unknown
							int hap=2;							
							if (scoreDiff > args.minScoreDifference) {
								allele = 1;
							}
							else if (scoreDiff < -args.minScoreDifference) {
								allele = 0;
							}
							else {
								allele = 2;
							}
							if (allele != 2) {

								chromSnps[curVar].GetPhase(allele, hap);
								if (args.verbosity > 0) {
									cerr << "allele: " << allele << " phase: " << hap << endl;
								}
								if (hap == 0) {
									h0++;
								}
								else if (hap == 1) {
									h1++;
								}
							}
							else {
								h2++;
							}
							
							if (args.verbosity > 0) {
								cerr << "pos: " << i << "\tscore: " << refScore << " " << altScore << " " << scoreDiff << " assigned: " << allele << "  nVars: " << endVar - startVar << " hap: " <<  hap  << " tpos: " << chromSnps[curVar].pos << endl;
							}
					}

					if (tAlnStr[i] != '-') {
						tPos ++;
					}
					if (qAlnStr[i] != '-') {
						qPos ++;
					}
				}
			}	
			
			bool haplotypeFound = false;
			if (args.assumeAutozygous and snvDb.size == 0) {
				haplotype = 2;
			}
			else {
				//				cout << h0 << " " << h1 << endl;
				if (args.verbosity > 0) {
					cerr << "h0: " << h0 << " h1: " << h1 << " minGenotyped " << args.minGenotyped << endl;
				}
				if (h0 >= args.minGenotyped or h1 >= args.minGenotyped) {
					if (h0 - args.minDifference > h1) {
						haplotype = 0;
						nH1+=1;
						bH1+= samRecord.seq.size();
						haplotypeFound = true;
						if (args.verbosity > 0) {
							cerr << "HAP0 "<< samRecord.title << endl;
						}
					}
					else if (h1 - args.minDifference > h0 ) {
						haplotype = 1;
						nH2+=1;
						bH2+= samRecord.seq.size();						
						haplotypeFound = true;
						if (args.verbosity > 0) {
							cerr << "HAP1 "<< samRecord.title << endl;
						}
					}
					else {
						haplotype = 2;
					}
				}
				if (haplotypeFound == false) {
					nUnk+=1;
					haplotype=2;
					bUnk+= samRecord.seq.size();											
				}
			}
			if (args.phaseStatsFileName != "") {
				phaseStatsOut << haplotype << "\t" << h0 << "\t" << h1 << "\t" << h2 << "\t" << refAlnLength << "\t"
											<< samRecord.chrom << ":" << refAlnStart << "-" << refAlnStart+refAlnLength << "\t"
											<< samRecord.title << "\t" << tPos + refAlnStart <<	endl;
			}
	
			if (splitHaplotypes) {
				if (haplotype == 0 || haplotype  == 1) {
					if (args.writeInterval == false) {
						haplotypeOut[haplotype] << samRecord.samLine << endl;
					}
					else {
						haplotypeOut[haplotype] << samRecord.chrom << "\t" << refAlnStart << "\t" << refAlnStart + refAlnLength << "\t" << samRecord.title << endl;
					}
				}
				else {
					if (args.unassigned != "") {
						if (args.writeInterval == false) {
							unassigned << samRecord.samLine << endl;
						}
						else {
							unassigned << samRecord.chrom << "\t" << refAlnStart << "\t" << refAlnStart + refAlnLength << "\t" << samRecord.title << endl;
						}
					}
					else {
						if (args.writeInterval == false) {
							int h;
							for (h=0;h<2;h++) { haplotypeOut[h] << samRecord.samLine << endl;	}
						}
						else {
							int h;
							for (h=0;h<2;h++) { haplotypeOut[h] << samRecord.chrom << "\t" << refAlnStart << "\t" << refAlnStart + refAlnLength << "\t" << samRecord.title << endl;}
						}
					}
				}
			}
			else {
				// add the haplotype tag to this file.
				if (args.writeInterval == false) {
					haplotypeOut[0] << samRecord.samLine << "\t" << args.phaseTag << ":Z:" << (int)haplotype
													<< "," << h0  << "," << h1 <<"," << h2  << endl;
				}
			}
		}
	}
	cerr << "partitionByPhasedSNVs ending with " << nH1 << "/" << bH1 << " " << nH2 << "/"
			 << bH2 << " " << nUnk << "/" << bUnk << endl;

	if (args.summaryFile != "") {
		ofstream summaryFile(args.summaryFile.c_str());
		summaryFile << "H1\t" << nH1 << "\tH2\t" << nH2
								<< "\tnUnk\t" << nUnk << "\tsites:\t"
								<< snvDb.size << endl;
	}

	
}

