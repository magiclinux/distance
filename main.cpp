#include "GraphUtil.h"
#include "Distance.h"
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cstring>

//#define TEST

// how to use the command line
//./distance -s 2 -n 10 -d -f edge_9.txt -k node_9.txt -w 0 -q newpool100.txt -R 3


static void usage() {
	cout << "\nUsage:\n"
		"	distance [-h] [-d] [-b] [-a alpha] [-t tc_limit] [-g n c r] [-n num_query] filename\n"
		"Description:\n"
		"	-h	Print the help message.\n"
		"	-g	Generate random graph (number of vertices, density, ratio).\n"
		"	-b	Perform shortest distance query by BFS.\n"
		"	-s	Build Spanning Tree Algorithm type(st-2 is default algorithm).\n"
		"		1:	st-1 (BFS Tree)\n"
		"		2:	st-2 (Maximal weight ST).\n"
		"	-a	Initial potential cost estimation alpha (default is 1.0).\n"
		"	-t	Initial TC_LIMIT for each iteration (default is 20,000,000).\n"
		"	-n	Set the total number of random queries. The default value is 100,000.\n"
		"	-d	debug mode (test distance query).\n"
		"	-f      filename	The graph file.\n"
		"	-k	keyword file.\n"
		<< endl;
}

static void keepResult(char* resultFileName, const char* version, char* filename, double alpha, 
		double ra, double ra2, int iteration, long tc_limit, float lt, float qt, int isize) {
	ofstream out(resultFileName, ios_base::out|ios_base::app);
	out << version << "_a" << alpha << "_ra" << ra << "_rb" << ra2 << "_tc" << tc_limit << "_iter" << iteration << "||" 
		<< filename << "\t" << lt << "\t" << qt << "\t" << isize << endl;
	out.close();
}

#ifdef TEST
int main(int argc, char* argv[]) {
	if (argc == 1) {
		usage();
		return 1;
	}
	
	float time = 4.6479e+06;
	cout << Util::formatTime(time) << endl;
	
	char* filename = argv[1];
	ifstream infile(filename);
	if (!infile) {
		cout << "Error: Cannot open " << filename << endl;
		return -1;
	}
	Graph g(infile);
	cout << "#vertex size:" << g.num_vertices() << "\t#edges size:" << g.num_edges() << endl;	
	int gsize = g.num_vertices();
	TDVec sptree;
	vector<int> rank;
	SparseVec treedist;
	GraphUtil::buildMaxWST(g, sptree, treedist, rank);
	// for test
	Util::printTDVec(sptree);
	cout << "All pairs shortest paths in spanning tree" << endl;
	Util::printSparseVec(treedist);
	for (int i = 0; i < rank.size(); i++)
		cout << rank[i] << "\t";
	cout << endl;
	
	int vid = 0;
	ReducedGraph treeedges;
	map<int,set<int> > right;
	map<int,set<int> >::iterator mit;
	set<int>::iterator sit;
	vector<int> left;
	int lsize = GraphUtil::genCandLabels(g, vid, sptree, treeedges, left, right);
	cout << "lsize = " << lsize << endl;
	cout << "-------------------TreeEdges-------------------" << endl;
	Util::printReducedGraph(treeedges);
	cout << "-------------------Left Label-------------------" << endl;
	for (int i = 0; i < left.size(); i++)
		cout << left[i] << " ";
	cout << endl;
	cout << "-------------------Right Label-------------------" << endl;
	for (mit = right.begin(); mit != right.end(); mit++) {
		cout << "v " << mit->first << ":";
		for (sit = mit->second.begin(); sit != mit->second.end(); sit++)
			cout << *sit << " ";
		cout << endl;
	}
}
#else
int main(int argc, char* argv[]) {
	if (argc == 1) {
		usage();
		return 1;
	}

	int i = 1;
	int st_type = 2;
	int query_num = 100000;
	long tc_limit = 20000000; 
	bool debug = false;
	bool grandom = false;
	int wlabel= 0;
	bool bfs = false;
	int vsize = 0, iter = 1;
	double density = 0, ratio = 1.0, alpha = 1.0;
	char *filename;
	char *keyfilename;
	char *queryfilename;
	int R;
	while (i < argc) {
		if (strcmp("-h", argv[i]) == 0) {
			usage();
			return 1;
		}
		if (strcmp("-g", argv[i]) == 0) {
			i++;
			grandom = true;
			vsize = atoi(argv[i++]);
			density = atof(argv[i++]);
			ratio = atof(argv[i++]);
		}
		if (strcmp("-d", argv[i]) == 0) {
			i++;
			debug = true;
		}	
		if (strcmp("-b", argv[i]) == 0) {
			i++;
			bfs = true;
		}			
		if (strcmp("-n", argv[i]) == 0) {
			i++;
			query_num = atoi(argv[i++]);
		}
		else if (strcmp("-s", argv[i]) == 0) {
			i++;
			st_type = atoi(argv[i++]);
		}
		else if (strcmp("-t", argv[i]) == 0) {
			i++;
			tc_limit = atoi(argv[i++]);
		}
		else if (strcmp("-a", argv[i]) == 0) {
			i++;
			alpha = atof(argv[i++]);
		}
		else if (strcmp("-f", argv[i]) == 0){
			i++;
			filename = argv[i++];
		}
		else if (strcmp("-k", argv[i]) == 0){
			i++;
			keyfilename = argv[i++];
		}
		else if (strcmp("-w", argv[i]) == 0){
			i++;
			wlabel = atoi(argv[i++]);	
			}
		else if (strcmp("-q", argv[i]) == 0){
			i++;
			queryfilename = argv[i++];
		}
		else if (strcmp("-R", argv[i]) == 0){
			i++;
			R = atoi(argv[i++]);
		}
	}
	
	if (grandom) {
		GraphUtil::genRandomRatioDAG(vsize, density, filename, ratio);
		exit(0);
	}
	
	ifstream infile(filename);
	/*
	if (!infile) {
		cout << "Error: Cannot open " << filename << endl;
		return -1;
	}
	*/
	cout << "INPUT: query_num=" << query_num << " st_type=" << st_type << " alpha=" << alpha << " tc_limit=" << tc_limit 
			<< " filename=" << filename << endl;
	
	Graph g((string)filename);
	//Graph g(infile);
	cout << "#vertex size:" << g.num_vertices() << "\t#edges size:" << g.num_edges() << endl;

	int s, t, dist;
	int left = 0;
	int gsize = g.num_vertices();

	// add by zhao
	g.readKeyword(keyfilename);
	//
	
	bool r;
	struct timeval after_time, before_time;
	float labeling_time = 0, query_time;

	// generate queries
	srand48(time(NULL));
	cout << "generating queries..." << endl;
	vector<int> src;
	vector<int> trg;
	vector<int>::iterator sit, tit;
	while (left < query_num) {
		s = lrand48() % gsize;
		t = lrand48() % gsize;
		src.push_back(s);
		trg.push_back(t);
		++left;
	}
	
	// perform BFS distance queries
	if (bfs) {
		gettimeofday(&before_time, NULL);
		for (sit = src.begin(), tit = trg.begin();
			sit != src.end(); ++sit, ++tit) {
			dist = GraphUtil::BFSDist(g,*sit,*tit);
		}
		gettimeofday(&after_time, NULL);
		query_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
			(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
		cout << "#total query running time:" << query_time << " (ms)" << endl;	
		keepResult("./result.txt", "BFS", filename, 0, 1, 0, 0, 0, 0, query_time, 0);
		exit(0);
	}
	
	string str = "TreeHop";
	if (st_type==1)
		str += "_ST1";
	else
		str += "_ST2";
		
	Distance dq(g);
	cout << "starting creating labels..." << endl;
	gettimeofday(&before_time, NULL);
	if(wlabel)
	{
		dq.createLabels(st_type,alpha,tc_limit);
		
	}
	else
	{
		dq.readAll();
	}
	gettimeofday(&after_time, NULL);
	labeling_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
		(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
	cout << "#labeling time:" << labeling_time << " (ms)\t Time " << Util::formatTime(labeling_time) << endl;


	cout << "process the kewang's algorithm"<<endl;
	
	FILE *fp = fopen(queryfilename,"r");
	int n;
	fscanf(fp,"%d",&n);
	for(int i=0; i<n; i++)
	{
		vector<string> v;
		int m; char buff[1000];
		fscanf(fp,"%d",&m);
		for(int j=0; j<m; j++)
		{
			fscanf(fp,"%s",buff);
			v.push_back((string)(buff));
		}
		gettimeofday(&before_time, NULL);
		dq.run(v,R);
		gettimeofday(&after_time, NULL);
		query_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
			(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
		cout << "#total query running time:" << query_time << " (ms)" << endl;
	}
	// process queries
	//cout << "process distance queries..." << endl;
/*
	gettimeofday(&before_time, NULL);
	
	if (debug) {
		cout << "queries test..." << endl;
		gettimeofday(&before_time, NULL);
		for (sit = src.begin(), tit = trg.begin(); sit != src.end(); ++sit, ++tit) {
			r = dq.test_distance(*sit, *tit);
		}
	}
	else {
		cout << "process queries..." << endl;
		gettimeofday(&before_time, NULL);
		for (sit = src.begin(), tit = trg.begin(); sit != src.end(); ++sit, ++tit) {
			dist = dq.distance(*sit, *tit);
		}
	}
	gettimeofday(&after_time, NULL);
	query_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
		(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
	cout << "#total query running time:" << query_time << " (ms)" << endl;
*/ //we don't need this part
	int labelsize, iteration;
	double ra = dq.stat_alpha();
	double ra2 = dq.stat_alpha2();
	labelsize = dq.label_size();
	cout << "#Label size = " << labelsize << endl;
	cout << "#Total iteration = " << dq.num_iter() << endl;
	keepResult("./result.txt", str.c_str(), filename, alpha, ra, ra2, iteration, tc_limit,
			labeling_time, query_time, labelsize);
}
#endif
