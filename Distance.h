#ifndef _DISTANCE_H
#define _DISTANCE_H

#include "GraphUtil.h"
#include "BiCover.h"

//#define DEBUG

class Distance {
	private:
		Graph& g;
		SparseVec dist;
		SparseVec treedist;
		vector<map<int,int> > l_in, l_out;

		// for test
		int iteration;
		double real_alpha, real_alpha2;
		struct timeval after_time, before_time;
		float run_time;	
		int *linearTree;

	public:
		Distance(Graph& graph);
		~Distance();
		void createLabels(int,double,long);
		void updateLabels(const vector<IDPair>& left, const vector<IDPair>& right);
		void addCandLabels(const vector<int>& rank, int start_id, int& end_id, 
				ReducedGraph& sptree, BiCover& bc, long tc_limit);
		void reducedGraphToLinear(const ReducedGraph& rg,int size);		
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
