#ifndef GENOTYPED_READ_H_
#define GENOTYPED_READ_H_



class GenotypedRead {
public:
	string chrom;
	int alnStart, alnEnd;
	vector<int> pos;
	vector<char> genotype;// 0, 1, or 2
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

bool ReadsOverlap(GenotypedRead &lhs, GenotypedRead &rhs) {
	return (lhs.chrom == rhs.chrom and 
					(lhs.alnStart < rhs.alnStart and lhs.alnEnd > rhs.alnStart) or
					(lhs.alnStart < rhs.alnEnd and   lhs.alnEnd > rhs.alnEnd) or
					(lhs.alnStart > rhs.alnStart and lhs.alnEnd < rhs.alnEnd));
}


#endif
