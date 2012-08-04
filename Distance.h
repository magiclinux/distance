#ifndef _DISTANCE_H
#define _DISTANCE_H

#include "GraphUtil.h"
#include "BiCover.h"

//#define DEBUG

struct y_element{
	int u;
	int d;
};

struct in_element{
	int x;
	int d;
};

struct out_element{
	int y;
	int d;
};

typedef y_element x_element;
struct y_invert_list{
	vector<y_element> l;
	map<int,int> m;
};
typedef y_invert_list x_invert_list;
struct entry{
	int u;
	int d;
	
};


class Distance {
	private:
		Graph& g;
		SparseVec dist;
		SparseVec treedist;
		vector<map<int,int> > l_in, l_out;

		map<int,map<int,int> > global_pair;// I can't remember why do we need this value.
		map<int,vector<entry> > global_invert;// also can't remember
		// for test
		int iteration;
		double real_alpha, real_alpha2;
		struct timeval after_time, before_time;
		float run_time;	
		int *linearTree;
		string labelfile;
		string sptreefile;

	public:
		Distance(Graph& graph);
		~Distance();
		void createLabels(int,double,long);
		void updateLabels(const vector<IDPair>& left, const vector<IDPair>& right);
		void addCandLabels(const vector<int>& rank, int start_id, int& end_id, 
				ReducedGraph& sptree, BiCover& bc, long tc_limit);
		void reducedGraphToLinear(const ReducedGraph& rg,int size);//add by zhao
		map<int,y_invert_list> constructY(string key,int D);// add by zhao
		void run(vector<string> keys,int D);		
		vector<in_element> getIn(int u, int D);
		vector<out_element> getOut(int u, int D);
		int nearestAncestor(int y, map<int,int> refs);
		bool check(vector<int> a, int v, int D);
		vector<vector<int> > step2(string key, int D, map<int,int> nodeAnswer, vector<vector<int> > answer );
		map<int,int> answerToNodeAnswer(vector<vector<int> > answer);
		
		void printNodeAnswer(map<int,int> nodeAnswer);
		void printAnswer(vector<vector<int> > answer);

		void printY(map<int,y_invert_list > Y);
		void printX(map<int,x_invert_list > X);

		void printGlobalPair();
		void constructGlobalInvert();
		void printGlobalInvert();
		
		void readAll();
		void writeAll();
		
		void writeLabels();
		void readLabels();
		
		void writeSptree(const ReducedGraph& rg);
		void readSptree();
		
		int distance(int, int);
		bool test_distance(int, int);
		int label_size();
		int num_iter();
		double stat_alpha();
		double stat_alpha2();
		
		// for test
		void showLabels();
};

#endif
