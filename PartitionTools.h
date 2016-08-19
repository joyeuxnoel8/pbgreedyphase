#ifndef PARTITIONING_TOOLS_H_
#define PARTITIONING_TOOLS_H_

#include <string>
using namespace std;

bool GetUngappedPrefix(string s, int p, int w, string &res) {
	// prefix is before p
	int i = p-1;
	int l = 0;
	while (i > 0 and l < w) {
		if (s[i] != '-') { l+=1;};
		i-=1;
	}
	if (l == w) {
		for (; i < p; i++) {
			if (s[i] != '-') { res.push_back(s[i]);}
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
	int i = p+1;
	int l = 0;
	while (i < s.size() and l < w) {
		if (s[i] != '-') { l+=1;};
		i+=1;
	}
	
	if (l == w) {
		for (; p < i; p++) {
			if (s[p] != '-') { res.push_back(s[p]);}
		}
		return true;
	}
	else {
		res = "";
		return false;
	}
}


#endif
