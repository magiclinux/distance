#ifndef UTIL_H_
#define UTIL_H_

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <deque>
#include <algorithm>
#include <utility>
#include <cmath>
#include <climits>
#include <ext/hash_map>
#include <sys/time.h>
#include <sys/resource.h>

using namespace std;
using namespace __gnu_cxx;

typedef map<int,vector<int> > 	ReducedGraph;
typedef map<int,map<int,int> >  SparseMatrix;
typedef vector<map<int,int> > 	SparseVec;
typedef vector<vector<int> > 	TDVec;
typedef vector<set<int> >		VecSet;
typedef unsigned int LINT; // it means the representation of largest number is 65535

template<class T>
inline string to_string(const T& t) {
	stringstream ss;
	ss << t;
	return ss.str();
}

class Util {
	public:
		static void printTDVec(const TDVec& tv);
		static void printReducedGraph(const ReducedGraph& rg);
		static void printSparseVec(const SparseVec& sv);
		static void printVecSet(const vector<set<int> >& vs);
		static void printPairVec(const vector<map<int,bool> >& vm);
		static void printMap(const map<LINT,int>& data);
		static string formatTime(float mstime);
		static void process_mem_usage(double& vm_usage, double& resident_set);
};

#endif
