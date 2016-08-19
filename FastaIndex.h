#ifndef FAI_H_
#define FAI_H_
#include <map>
#include <string>
#include <vector>
using namespace std;


class FastaIndex {
 public:
	map<string, vector<long> > fai;
	ifstream fastaIn;
	
	void Initialize(string name) {
		fastaIn.open(name.c_str());
		string faiName = name + ".fai";
		ifstream faiIn(faiName.c_str());
		string contig;
		long length, start, lineSeq, lineLen;
		while ((faiIn >> contig >> length >> start >> lineSeq >> lineLen)) {
			fai[contig] = vector<long>(4);
			fai[contig][0] = length;
			fai[contig][1] = start;
			fai[contig][2] = lineSeq;
			fai[contig][3] = lineLen;
		}
	}
	void Print() {
		return;
		map<string, vector<long> >::iterator faiIt;
		for (faiIt = fai.begin(); faiIt != fai.end(); ++faiIt) {
			int i;
			cout << (*faiIt).first << " ";
			for (i = 0; i < 5; i++) {
				cout << (*faiIt).second[i] << " ";
			}
			cout << endl;
		}
	}

	void GetSeq(string &seq, string &contig, int start=0, int end=0) {
	  if (fai.find(contig) == fai.end()) {
		start=0; 
		end = 0;
		seq="";
		return;
	  }
		this->Print();
		vector<long> &v = fai[contig];
		long lineLen = fai[contig][2];
		long fileLineLength = fai[contig][3];
		long chrStart = fai[contig][1];
		long seqLength;
		if (start == end and start == 0) {
			start = 0;
			end = fai[contig][0];
		}
		long startLine = start / lineLen;
		long endLine   = end   / lineLen;
		long startLinePos = start % lineLen;
		long endLinePos = end % lineLen;
		

    long startFilePos = chrStart + startLine * fileLineLength + startLinePos;
		long endFilePos   = chrStart + endLine * fileLineLength  + endLinePos;

		long fileReadLength = endFilePos - startFilePos;

		fastaIn.seekg(startFilePos);
		string fileBuf;
		fileBuf.resize(endFilePos - startFilePos);

		fastaIn.get(&fileBuf[0], endFilePos - startFilePos, '\0');
		seq.resize(end - start);
		int i,j;
		j = 0;
		for (i = 0; i < fileBuf.size(); i++) {
			if (fileBuf[i] != '\n' and fileBuf[i] != '\0') {
				assert(j < seq.size());
				seq[j] = fileBuf[i];
				j++;
			}
		}
	}
};
	


#endif
