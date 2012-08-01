#ifndef _GRAPH_UTIL_H_
#define _GRAPH_UTIL_H_

#include "Graph.h"

#define ALLPATHS

class GraphUtil {
	public:
		// for scalable version
		static void buildMaxWST(Graph& g, ReducedGraph& sptree, vector<int>& rank);
		static void buildSPTree(Graph& g, ReducedGraph& sptree, vector<int>& rank);
		
		// directly generate candidate labels for SCP on each BFS tree
		static int genCandLabels(Graph& g, int vid, ReducedGraph& sptree, 
				map<int,set<int> >& tpairs, vector<int>& left, map<int,set<int> >& right);
		
		// start_id and end_ind are indices of rank vector
		static void computePartialSD(Graph& g, vector<int>& rank, int start_ind, 
				int& end_ind, SparseVec& pdist, long tc_limit); 
		static void computeShortestDistanceOnTree(ReducedGraph& sptree, SparseVec& dist, int gsize);
	
		// for efficient version
		static void computeShortestDistance(Graph& g, SparseVec& dist);
		static int  BFSTree(Graph& g, int v, ReducedGraph& tree);
		static void PVisit(ReducedGraph& tree, int v, int& id, int& cn, hash_map<int,int>& order, hash_map<int,int>& min);
		static int  buildSpanningTree1(Graph& g, SparseVec& treedist, SparseVec& dist, ReducedGraph& sptree);
		static void treestat(ReducedGraph& tree, map<int,pair<int,int> >& pred_succ, int root, int gsize);
		static void treetc(ReducedGraph& tree, map<int,set<int> >& tpairs);
		static void treetcvisit(ReducedGraph& tree, int root, map<int,set<int> >& tpairs, map<int,bool>& mark);
		static void treetc1(ReducedGraph& tree, int root, map<int,set<int> >& tpairs);
		static void distPairs(ReducedGraph& tree, int root, SparseVec& treedist, SparseVec& dist);
		static bool checkPath(Graph& g, int sid, int tid);
		static int  BFSDist(Graph& g, int sid, int tid);
		static int  buildSpanningTree2(Graph& g, SparseVec& treedist, SparseVec& dist, ReducedGraph& sptree);
		
		// for test
		static void genRandomGraph(int n, double c, char* filename);
		static void genRandomDAG(int n, double c, char* filename);
		static void genRandomRatioDAG(int n, double c, char* filename, double ratio);

};

#endif
