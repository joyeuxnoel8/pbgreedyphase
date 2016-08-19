#ifndef SNVDB_H_
#define SNVDB_H_

class SNV {
public:
	int pos;
	char nuc, ref;
	int h1, h2;
	int nalt, nref;
	SNV(int p) {
		pos = p;
		nuc = '\0';
		ref = '\0';
		h1 = h2 = 0;
		nalt = nref = 0;
	}

	SNV(int p, char n, char r, int _h1, int _h2) {
		pos = p;
		nuc = n;
		ref = r;
		h1  = _h1;
		h2  = _h2;
		nalt = nref = 0;
	}

	void Print() {
		cout << pos << "\t" << ref << "\t" << nuc << "\t" << nalt << "\t" << nref << endl;
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
		nalt = rhs.nalt;
		nref = rhs.nref;
		return *this;
	}
	
	bool GetPhase(int allele, int &phase) {
		if (allele == h1) { phase = 0; return true; }
		if (allele == h2) { phase = 1; return true; }
		phase = 2;
		return false;
	}
};



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

	
	void AddSNV(string chrom, int pos, // location
							char nuc,  // alt
							char ref,  // ref
							int h1=0, int h2=0) {
		
		if (db.find(chrom) == db.end()) {
			db[chrom] = SNVs();
		}

		db[chrom].push_back(SNV(pos, nuc, ref, h1, h2));
		size++;
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

	void Print() {
		ChromDB::iterator chromIt, chromEnd;
		for (chromIt = db.begin(); chromIt != db.end(); ++chromIt) {
			int i;
			for (i = 0; i < chromIt->second.size(); i++) {
				chromIt->second[i].Print();
			}
		}
	}
};
#endif
