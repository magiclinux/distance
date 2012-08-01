#include "Util.h"

void Util::printTDVec(const TDVec& tv) {
	vector<int>::const_iterator vit;
	cout << "TDVec" << endl;
	for (int i = 0; i < tv.size(); i++) {
		cout << "v " << i << ":";
		for (vit = tv[i].begin(); vit != tv[i].end(); vit++)
			cout << *vit << " ";
		cout << endl;
	}
}

void Util::printSparseVec(const SparseVec& sv) {
	SparseVec::const_iterator sit;
	map<int,int>::const_iterator mit;
	cout << "SparseVec" << endl;
	for (int i = 0; i < sv.size(); i++) {
		cout << "from " << i << ": ";
		for (mit = sv[i].begin(); mit != sv[i].end(); mit++) 
			cout << mit->first << "[" << mit->second << "] ";
		cout << endl;
	}
}

void Util::printReducedGraph(const ReducedGraph& rg) {
	ReducedGraph::const_iterator rit;
	vector<int>::const_iterator vit;
	cout << "ReducedGraph " << endl;
	for (rit = rg.begin(); rit != rg.end(); rit++) {
		cout << "vertex " << rit->first << ": ";
		for (vit = rit->second.begin(); vit != rit->second.end(); vit++)
			cout << *vit << " ";
		cout << " #" << endl;
	}
	cout << endl;
}

void Util::printVecSet(const vector<set<int> >& vs) {
	set<int>::const_iterator sit;
	for (int i = 0; i < vs.size(); i++) {
		cout << i << " : ";
		for (sit = vs[i].begin(); sit != vs[i].end(); sit++) 
			cout << *sit << " ";
		cout << endl;
	}
	cout << endl;
}

void Util::printPairVec(const vector<map<int,bool> >& vm) {
	map<int,bool>::const_iterator mit;
	for (int i = 0; i < vm.size(); i++) {
		cout << "v " << i << " : ";
		for (mit = vm[i].begin(); mit != vm[i].end(); mit++)
			cout << "[" << mit->first << "," << mit->second << "] ";
		cout << endl;
	}
	cout << endl;
}

void Util::printMap(const map<LINT,int>& data) {
	map<LINT,int>::const_iterator mit;
	for (mit = data.begin(); mit != data.end(); mit++) {
		cout << "{" << mit->first << "," << mit->second << "} ";
	}
	cout << endl;
}

string Util::formatTime(float mstime) {
	int mtime = (int)mstime;
	string timestr = "";
	int unit[] = {60, 60, 60, 1000};
	float quality[] = {0, 0, 0, 0};
	string ustr[] = {"Hours", "Mins", "Secs", "Mss"};
	
	int count = 3;
	while (count>=0) {
		quality[count] = mtime%unit[count];
		mtime = (int)(mtime/unit[count]);
		count--;
	}
	
	timestr = to_string(quality[0])+":"+to_string(quality[1])+":"+to_string(quality[2])
			+" (" +to_string(quality[3])+"ms)";
	return timestr;
}

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
void Util::process_mem_usage(double& vm_usage, double& resident_set) {
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

