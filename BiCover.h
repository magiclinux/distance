#ifndef _BICOVER_H
#define _BICOVER_H

#include "GraphUtil.h"

#define MFACTOR 10

//#define DEBUG

#define STAT_ALPHA

typedef pair<LINT,LINT>  LabelPair;
typedef pair<int,LINT>   IntPair;
typedef pair<int,int>	 IDPair;

struct strictPairComp{
	bool operator() (const IntPair& l1, const IntPair& l2) const {
		if (l1.first<l2.first || (l1.first==l2.first && l1.second<l2.second))
			return true;
		return false;
	}
};

typedef map<IntPair,int,strictPairComp>			PairMap;
typedef set<IntPair,strictPairComp> 			PairSet;
typedef vector<map<int,bool> >					PairVec;
typedef vector<map<int,vector<LabelPair> > >	VecMap;

class BiCover{
	private:
		map<LINT,list<LINT> > leftlist, rightlist;
		VecMap spmap;
		PairVec covered;
		map<LINT,int> leftcover, rightcover;
		PairSet leftps, rightps;
		
		map<LINT,bool> leftselected, rightselected;
		set<LINT>	leftsingle, rightsingle;
		int gsize, lpsize, rpsize, num_covered, tlsize;
		double alpha, real_alpha2;
		
		// for test
		struct timeval after_time, before_time;
		float run_time;			

		// for compute alpha
		#ifdef STAT_ALPHA
		map<LINT,int> lsv, ltv, rsv, rtv;
		#endif
		
	public:
		BiCover(int gs, double init_alpha);
		~BiCover();
		void addEdge(int uid, int xid, int yid, int vid);
		void scp(int cycleid, int end_ind);
		void constructMapping(int cycleid);
		pair<LINT,bool> selectObj(int sideno);
		void updateState(LINT selectid, bool leftmark);
		void update();
		int update(vector<IDPair>& left, vector<IDPair>& right);
		void clear();
		vector<IDPair> leftPairs();
		vector<IDPair> rightPairs();		
		void printBipartiteGraph();
		void printPairSet(const PairSet& ps);
		double stat_alpha();
		double stat_alpha2();
};

#endif
