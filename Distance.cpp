#include "Distance.h"

/*
*************************************************************************
* Implementation of class Distance
*************************************************************************
*/

Distance::Distance(Graph& graph):g(graph) {
	l_in = vector<map<int,int> >(g.num_vertices(),map<int,int>());
	l_out = vector<map<int,int> >(g.num_vertices(),map<int,int>());
	real_alpha2 = 0;
}

Distance::~Distance() {
	delete [] linearTree;
}

void Distance::createLabels(int sttype, double alpha, long tc_limit) {
	struct timeval stime, etime;
	float  rtime;

	// build spanning tree
	ReducedGraph sptree;
	vector<int> rank;
	int gsize = g.num_vertices();
	gettimeofday(&before_time, NULL);
	if (sttype==1) {
		GraphUtil::buildSPTree(g, sptree, rank);
	}
	else {
		GraphUtil::buildMaxWST(g, sptree, rank);
	}
	gettimeofday(&after_time, NULL);
	run_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
		(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
	cout << "#build spanning tree computation time:" << run_time << " (ms)" << endl;	
	
	// main loop
	int start_ind = 0, end_ind = 0;
	int rm_labelsize = 0;
	long inc_covered = 0;
	vector<IDPair> left, right;
	BiCover bc = BiCover(gsize,alpha);

	iteration = 0;
	do {
		// clear temp labels
		vector<IDPair>().swap(left);
		vector<IDPair>().swap(right);
	
		gettimeofday(&stime, NULL);
	
		// add candidate labels
		gettimeofday(&before_time, NULL);
		addCandLabels(rank, start_ind, end_ind, sptree, bc, tc_limit);
		gettimeofday(&after_time, NULL);
		run_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
			(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
		cout << "\tCycleID " << iteration << "\taddCandLabels computation time:" << run_time << " (ms)" << endl;	

		// for test
		#ifdef DEBUG
		cout << "end_ind=" << end_ind << endl;
		bc.printBipartiteGraph();
		#endif

		// init BiCover
		gettimeofday(&before_time, NULL);
		bc.constructMapping(iteration);
		gettimeofday(&after_time, NULL);
		run_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
			(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
		cout << "\tCycleID " << iteration << "\tconstructMapping computation time:" << run_time << " (ms)" << endl;	
		
		// perform greedy set cover
		gettimeofday(&before_time, NULL);
		bc.scp(iteration,end_ind);
		gettimeofday(&after_time, NULL);
		run_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
			(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
		cout << "\tCycleID " << iteration << "\tSCP computation time:" << run_time << " (ms)" << endl;		

		// update BiCover after each scp operation
		gettimeofday(&before_time, NULL);
		rm_labelsize = bc.update(left, right);
		gettimeofday(&after_time, NULL);
		run_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
			(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
		cout << "\tCycleID " << iteration << "\tupdateSCP computation time:" << run_time << " (ms)" << endl;

		// update treehop labels
		gettimeofday(&before_time, NULL);
		updateLabels(left, right);
		gettimeofday(&after_time, NULL);
		run_time = (after_time.tv_sec - before_time.tv_sec)*1000.0 + 
			(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
		cout << "\tCycleID " << iteration << "\tupdateLabels computation time:" << run_time << " (ms)" << endl;

		// update tc_limit based on covering progress
		inc_covered = left.size()+right.size();
		/*
		if (inc_covered>0)
			tc_limit = inc_covered;
		else
			tc_limit = tc_limit/3;
		*/
		
		// update tc_limit based on left label size
		if (tc_limit == rm_labelsize)
			tc_limit = tc_limit/4;
		else if (tc_limit > rm_labelsize)
			tc_limit = tc_limit-rm_labelsize;
			
		cout << "\t# current_labels=" << inc_covered << "\ttc_limit=" << tc_limit << "\tleft_labelsize=" 
				<< rm_labelsize << "\tstart_ind=" << start_ind << "\tend_ind=" << end_ind << endl;
		#ifdef DEBUG
		showLabels();
		#endif
							
		// for test
/* 		if (iteration == 1)
		exit(0);	 */	
/* 		char ch;
		cin >> ch; */
		
		gettimeofday(&etime, NULL);
		rtime = (etime.tv_sec - stime.tv_sec)*1000.0 + (etime.tv_usec - stime.tv_usec)*1.0/1000.0;
		cout << "##CycleID " << iteration << " total computation time:" << rtime << " (ms)\tTime " 
				<< Util::formatTime(rtime) << "\tTotal_size=" << label_size() << endl;	
						
		start_ind = end_ind;
		iteration++;
	} while (start_ind < gsize);
	
	#ifdef STAT_ALPHA
	real_alpha = bc.stat_alpha();
	cout << "# real alpha = " << real_alpha << endl;
	real_alpha2 = bc.stat_alpha2();
	cout << "# real_alpha2 = " << real_alpha2 << endl;
	#endif
	
	// clear all data structures
	bc.clear();
	vector<int>().swap(rank);
	vector<IDPair>().swap(left);
	vector<IDPair>().swap(right);	
	
	//#ifdef DEBUG
	Util::printReducedGraph(sptree);
	//#endif
	reducedGraphToLinear(sptree,g.num_vertices());// add by zhao
	
	// compute all pairs shortest paths for sptree
	GraphUtil::computeShortestDistanceOnTree(sptree,treedist,gsize);
	ReducedGraph().swap(sptree);
	
	// for test
	//#ifdef DEBUG
	showLabels();
	Util::printSparseVec(treedist);
	//#endif
}

void Distance::reducedGraphToLinear(const ReducedGraph& rg,int size) {// add by zhao
	ReducedGraph::const_iterator rit;
	vector<int>::const_iterator vit;
	cout << "ReducedGraph " << endl;
	linearTree = new int[size];
	for(int i=0; i<size; i++)
		linearTree[i]=-1;
	for (rit = rg.begin(); rit != rg.end(); rit++) {
		//cout << "vertex " << rit->first << ": ";
		for (vit = rit->second.begin(); vit != rit->second.end(); vit++)
			linearTree[*vit] = rit->first;
			//cout << *vit << " ";
		//cout << " #" << endl;
	}
	cout << "Look at the LinearTree:"<<endl;
	for(int i=0; i<size; i++)
		cout<< "["<<i<<"]="<<linearTree[i]<<endl;
	cout<<endl;
	//cout << endl;
}

void Distance::addCandLabels(const vector<int>& rank, int start_id, int& end_id, 
		ReducedGraph& sptree, BiCover& bc, long tc_limit) {
	map<int,set<int> > tpairs;
	map<int,set<int> > right;
	vector<int> left;
	
	ReducedGraph::iterator rit;
	vector<int>::iterator  vit;
	vector<int>::const_iterator svit;
	set<int>::iterator     sit;
	set<int>::const_iterator scit;
	map<int,set<int> >::iterator mit;
	
	long inc_count = 0, lsize;
	int j = start_id;
	int vid = rank[j];
	while (true) {
		// generate candidate labels for SCP
		// tpairs.clear();
		map<int,set<int> >().swap(tpairs);
		// left.clear();
		vector<int>().swap(left);
		// right.clear();
		map<int,set<int> >().swap(right);
	//	cout << "start genCandLabels" << endl;
		lsize = GraphUtil::genCandLabels(g, vid, sptree, tpairs, left, right);
	//	cout << "complete genCandLabels " << endl;
		
		#ifdef DEBUG
		cout << "vid=" << vid << "\tlsize = " << lsize << endl;
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
		#endif
		
		// add new labels into BiCover
		for (vit = left.begin(); vit != left.end(); vit++) {
			// hop is one vertex (special case of tree hop)
			for (sit = right[*vit].begin(); sit != right[*vit].end(); sit++) {
				if (vid != *sit)
					bc.addEdge(vid,*vit,*vit,*sit);
			}
			// regular case
			for (scit = tpairs[*vit].begin(); scit != tpairs[*vit].end(); scit++) {
				for (sit = right[*scit].begin(); sit != right[*scit].end(); sit++) {
					// for test
			//		cout << "add bicover edges (" << vid << "," << *vit << "," << *scit << "," << *sit << ")" << endl;
					bc.addEdge(vid,*vit,*scit,*sit);
				}
			}
		}
		
		inc_count += lsize;
		j++;
		if (inc_count>=tc_limit || j>=rank.size())
			break;
			
		vid = rank[j];
	}
	end_id = j;
}

// update distance labels
void Distance::updateLabels(const vector<IDPair>& lp, const vector<IDPair>& rp) {
	int sid, tid, sd;
	map<int, map<int,int> > dist;
	for (int i = 0; i < lp.size(); i++) {
		sid = lp[i].first;
		tid = lp[i].second;
		if (sid == tid)
			l_out[sid].insert(make_pair(tid,0));
		else if (l_out[sid][tid] == 0) {
			sd = GraphUtil::BFSDist(g,sid,tid);
			l_out[sid][tid] = sd;
//			dist[sid][tid] = sd;
		}
	}
	for (int i = 0; i < rp.size(); i++) {
		sid = rp[i].first;
		tid = rp[i].second;
		if (sid == tid)
			l_in[tid].insert(make_pair(sid,0));
		else if (l_in[tid][sid] == 0) {
			sd = GraphUtil::BFSDist(g,sid,tid);
			l_in[tid][sid] = sd;
//			dist[sid][tid] = sd;
		}
	}
}

int Distance::num_iter() {
	return iteration;
}

double Distance::stat_alpha() {
	return real_alpha;
}

double Distance::stat_alpha2() {
	return real_alpha2;
}

void Distance::showLabels() {
	map<int,int>::iterator mit;
	cout << "====================Distance Labels====================" << endl;
	cout << "OutLabels" << endl;
	for (int i = 0; i < l_out.size(); i++) {
		cout << "lout " << i << " : ";
		for (mit = l_out[i].begin(); mit != l_out[i].end(); mit++)
			cout << "[" << mit->first << "," << mit->second << "] ";
		cout << endl;
	}
	cout << "InLabels" << endl;
	for (int i = 0; i < l_in.size(); i++) {
		cout << "lin " << i << " : ";
		for (mit = l_in[i].begin(); mit != l_in[i].end(); mit++)
			cout << "[" << mit->first << "," << mit->second << "] ";
		cout << endl;
	}	
	cout << endl;
}

int Distance::label_size() {
	int size = 0;
	for (int i = 0; i < l_in.size(); i++)
		size += l_in[i].size();
	for (int i = 0; i < l_out.size(); i++)
		size += l_out[i].size();
	return size;
}

int Distance::distance(int src, int trg) {
	if (src == trg) return 0;
	
	int mindist = MAX_VAL;
	map<int,int>::iterator mit1, mit2;
	for (mit1 = l_out[src].begin(); mit1 != l_out[src].end(); mit1++) {
		if (mit1->first == trg) return mit1->second;
		for (mit2 = l_in[trg].begin(); mit2 != l_in[trg].end(); mit2++) {
			if (mit1->first == mit2->first) {
				if (mindist>(mit1->second+mit2->second)) {
					mindist = mit1->second+mit2->second;
				}
			}
			else if (treedist[mit1->first].find(mit2->first) != treedist[mit1->first].end()) {
				if (mindist>(mit1->second+mit2->second+treedist[mit1->first][mit2->first])) {
					mindist = mit1->second+mit2->second+treedist[mit1->first][mit2->first];
				}
			}
		}
	}
	if (mindist < MAX_VAL) return mindist;
	return -1; // unreachable
}	

bool Distance::test_distance(int src, int trg) {
	int ans = distance(src,trg);
	int sd = GraphUtil::BFSDist(g,src,trg);
	if (ans!=sd) {
		cout << "Wrong: [" << src << "] to [" << trg << "] distance = " << ans << endl;
		cout << "correct distance = " << sd << endl;
		exit(0);
	}
	return ans==sd;
}

