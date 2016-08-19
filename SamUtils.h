#ifndef SAM_UTIS_H_
#define SAM_UTIS_H_


class SAMHeader {
public:
	vector<string> lines;
	vector<string> chroms;
	vector<int> lengths;
};


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

void FixedNameToIndexMap(SAMHeader &header, vector<FASTASequence> &genome, map<string, int> &refIndex) {
	int i;
	for (i = 0; i < header.chroms.size(); i++) {
		refIndex[header.chroms[i]] = 0;
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

void MakeAlignStrings(SAMRecord &record, char *ref, int regionStart, int regionEnd, string &queryAln, string &refAln) {
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
		if (record.ops[i] == 'M' or record.ops[i] == 'X' or record.ops[i] == '=') {
			if (regionEnd != -1 and refPos - record.refPos > regionEnd - regionStart) {
				cerr << "ERROR! Alignment " << record.title << " extends past reference region." << endl;
				exit(0);
			}
			strncpy(&queryAln[queryAlnPos], &query[queryPos], record.lengths[i]);
			strncpy(&refAln[refAlnPos], &ref[refPos-regionStart], record.lengths[i]);
			queryAlnPos += record.lengths[i];
			queryPos += record.lengths[i];
			refAlnPos += record.lengths[i];
			refPos += record.lengths[i];

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
			if (regionEnd != -1 and refPos - record.refPos > regionEnd - regionStart) {
                                cerr << "ERROR! Alignment " << record.title << " extends past reference region." << endl;
                                exit(0);
                        }

			strncpy(&refAln[refAlnPos], &ref[refPos-regionStart], record.lengths[i]);
			refAlnPos += record.lengths[i];
			refPos += record.lengths[i];
			int j;
			for (j = 0; j < record.lengths[i]; j++, queryAlnPos++) {
				queryAln[queryAlnPos] = '-';
			}
		}
	}
}

#endif
