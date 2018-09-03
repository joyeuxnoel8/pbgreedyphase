#ifndef PARTITIONING_TOOLS_H_
#define PARTITIONING_TOOLS_H_

#include <string>
#include <vector>
#include <ctype.h>

using namespace std;

string JoinPrefixSuffix(string p, char c, string s ) {
	if (c != '-') {
		return p + c + s;
	}
	else {
		return p + s;
	}
}
																											
bool GetUngappedPrefix(string s, int p, int w, string &res) {
	// prefix is before p
	int i = p-1;
	int l = 0;
	res = "";
	while (i > 0 and l < w) {
		if (s[i] != '-') { l+=1;};
		i-=1;
	}
	if (l == w) {
		i+=1;
		for (; i <= p-1; i++) {
			if (s[i] != '-') { res += toupper(s[i]);}
		}
		// reverse the result
		for (i = 0; i < w/2; i++) {
			char tmp = res[i];
			res[i] = res[res.size()-i-1];
			res[res.size()-i-1] = tmp;
		}
		return true;
	}
	else {
		res = "";
		return false;
	}
}

bool GetUngappedSuffix(string s, int p, int w, string &res) {
	// prefix is before p
	res = "";
	int i = p+1;
	int l = 0;
	while (i < s.size() and l < w) {
		if (s[i] != '-') { l+=1;};
		i+=1;
	}
	
	if (l == w) {
		p+=1;
		for (; p < i; p++) {
			if (s[p] != '-') { res.push_back(toupper(s[p]));}
		}
		return true;
	}
	else {
		res = "";
		return false;
	}
}

int Index(int i, int j, int colLen) {
	return i*colLen + j;
}
int SWAlign(string q, string t, int match, int mismatch, int gap, string &qAlnStr, string &tAlnStr) {
	vector<int> score, path;
	int rowLen = q.size()+1;
	int colLen = t.size()+1;
	score.resize(rowLen*colLen, 0);
	path.resize(rowLen*colLen, 0);
	int i,j;
	for (i = 1; i < colLen;i++) {
		score[Index(i, 0, colLen)] = score[Index(i-1,0,colLen)]+gap;
		path[Index(i, 0, colLen)] = 1;
	}
	for (j = 1; j < rowLen; j++) {
		score[Index(0,j,colLen)] = score[Index(0,j-1,colLen)]+gap;
		path[Index(0,j,colLen)]  = 2;
	}
	for (i = 1; i < rowLen; i++) {
		for (j = 1; j < colLen; j++) {
			char qc = q[i-1];
			char tc = t[j-1];
			int ms, is, ds;
			if (qc == tc) {
				ms = score[Index(i-1,j-1,colLen)] + match;
			}
			else {
				ms = score[Index(i-1,j-1,colLen)] + mismatch;
			}
			is = score[Index(i-1,j,colLen)] + gap;
			ds = score[Index(i,j-1,colLen)] + gap;
			if (ms < is and ms < ds) {
				score[Index(i,j,colLen)] = ms;
				path[Index(i,j,colLen)] = 0;
			} else if (is < ms and is < ds) {
				score[Index(i,j,colLen)] = is;
				path[Index(i,j,colLen)] = 1;
			} else {
				score[Index(i,j,colLen)] = ds;
				path[Index(i,j,colLen)] = 2;
			}
			//			cout << setw(4) << score[Index(i,j,colLen)] << " ";
		}
		//		cout << endl;
	}
	i = colLen-1;
	j = colLen-1;
	qAlnStr = "";
	tAlnStr = "";
	while (i > 0 or j > 0) {
		if (path[Index(i,j,colLen)] == 0) {
			qAlnStr+=q[i-1];
			tAlnStr+=t[j-1];
			i-=1;
			j-=1;
		}
		else if (path[Index(i,j,colLen)] == 1) {
			qAlnStr+=q[i-1];
			tAlnStr+='-';
			i-=1;
		}
		else if (path[Index(i,j,colLen)] == 2) {
			qAlnStr+='-';
			tAlnStr+=t[j-1];
			j-=1;
		}
	}
	for (i = 0; i < qAlnStr.size()/2; i++) {
		char tmp = qAlnStr[i];
		qAlnStr[i] = qAlnStr[qAlnStr.size()-1-i];
		qAlnStr[qAlnStr.size()-1-i] = tmp;
	}
	for (i = 0; i < tAlnStr.size()/2; i++) {
		char tmp = tAlnStr[i];
		tAlnStr[i] = tAlnStr[tAlnStr.size()-1-i];
		tAlnStr[tAlnStr.size()-1-i] = tmp;
	}
	//	cout << tAlnStr << endl;
	//	cout << qAlnStr << endl;
	return score[Index(rowLen-1,colLen-1,colLen)];		
}

#endif
