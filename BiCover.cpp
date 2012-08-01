#include "BiCover.h"

BiCover::BiCover(int gs, double init_alpha): gsize(gs), alpha(init_alpha) {
	num_covered = 0;
	spmap = VecMap(gsize, map<int,vector<LabelPair> >());
	covered = PairVec(gsize, map<int,bool>());
	real_alpha2 = 0;
	tlsize = 0;
}

BiCover::~BiCover() {
}

void BiCover::clear() {
	// clear all data structures
	map<LINT,list<LINT> >().swap(leftlist);
	map<LINT,list<LINT> >().swap(rightlist);
	VecMap().swap(spmap);
	PairVec().swap(covered);
	map<LINT,int>().swap(leftcover);
	map<LINT,int>().swap(rightcover);
	PairSet().swap(leftps);
	PairSet().swap(rightps);
	map<LINT,bool>().swap(leftselected);
	map<LINT,bool>().swap(rightselected);
	set<LINT>().swap(leftsingle);
	set<LINT>().swap(rightsingle);	
}

void BiCover::addEdge(int uid, int xid, int yid, int vid) {
	LINT leftid = uid*gsize+xid;
	LINT rightid = vid*gsize+yid;
	leftlist[leftid].push_back(rightid);
	rightlist[rightid].push_back(leftid);
//	if (!leftselected[leftid])
	if (leftselected.find(leftid) == leftselected.end())
		leftcover[leftid] += (int)(MFACTOR*alpha);
	else {
		#ifdef STAT_ALPHA
		lsv[leftid]++;
		ltv[leftid]++;
		#endif
	}
//	if (!rightselected[rightid])
	if (rightselected.find(rightid) == rightselected.end())
		rightcover[rightid] += (int)(MFACTOR*alpha);
	else {
		#ifdef STAT_ALPHA
		rsv[rightid]++;
		rtv[rightid]++;
		#endif
	}
	// update covered
	covered[uid].insert(make_pair(vid,false));
	// update edge mapping for speedup
	spmap[uid][vid].push_back(make_pair(leftid,rightid));
}

// mapping leftcover rightcover to leftps and rightps
void BiCover::constructMapping(int cycleid) {
	// avoid redundant initilization
//	if (cycleid>0) return;
	
//	leftps.clear();
//	rightps.clear();
	PairSet().swap(leftps);
	PairSet().swap(rightps);
	map<LINT,int>::const_iterator mit;
	for (mit = leftcover.begin(); mit != leftcover.end(); mit++) {
		leftps.insert(IntPair(mit->second,(mit->first)));
	}
	for (mit = rightcover.begin(); mit != rightcover.end(); mit++) {
		rightps.insert(IntPair(mit->second,(mit->first)));
	}
	cout << "\t# of left label pairs:" << leftcover.size() << "\tright label pairs:" << rightcover.size() << endl;

	// for test
	#ifdef DEBUG
	cout << "----------------------leftps--------------------" << endl;
	printPairSet(leftps);
	cout << "----------------------rightps--------------------" << endl;
	printPairSet(rightps);
	#endif
}

void BiCover::update() {
	list<LINT>::iterator lit;
	map<int,bool>::iterator  miter;
	map<LINT,int>::iterator  mlit;
	map<LINT,bool>::iterator mbit;

	// update leftlist and rightlist according to leftcover and rightcover
	map<LINT,list<LINT> >  tmp_left, tmp_right;
	for (mlit = leftcover.begin(); mlit != leftcover.end(); mlit++) {
		for (lit = leftlist[mlit->first].begin(); lit != leftlist[mlit->first].end(); lit++) {
			if (!covered[mlit->first/gsize][*lit/gsize])
				tmp_left[mlit->first].push_back(*lit);
		}
	}
//	leftlist.clear();
	map<LINT,list<LINT> >().swap(leftlist);
	leftlist = tmp_left;
//	tmp_left.clear();
	map<LINT,list<LINT> >().swap(tmp_left);

	for (mlit = rightcover.begin(); mlit != rightcover.end(); mlit++) {
		for (lit = rightlist[mlit->first].begin(); lit != rightlist[mlit->first].end(); lit++) {
			if (!covered[*lit/gsize][mlit->first/gsize])
				tmp_right[mlit->first].push_back(*lit);
		}
	}
//	rightlist.clear();
	map<LINT,list<LINT> >().swap(rightlist);
	rightlist = tmp_right;
//	tmp_right.clear();	
	map<LINT,list<LINT> >().swap(tmp_right);

	// update covered
	for (int i = 0; i < covered.size(); i++) {
		if (covered[i].size() == 0) continue;
		for (miter = covered[i].begin(); miter != covered[i].end(); ) {
			if (miter->second) 
				covered[i].erase(miter++);
			else
				miter++;
		}
	}
		
	// clear spmap
	if (leftlist.size()+rightlist.size() == 0) {
		for (int i = 0; i < spmap.size(); i++)
			map<int,vector<LabelPair> >().swap(spmap[i]);
			
		// clear leftcover right cover leftps rightps
		map<LINT,int>().swap(leftcover);
		map<LINT,int>().swap(rightcover);
		PairSet().swap(leftps);
		PairSet().swap(rightps);
		
		// clear leftselected
		for (mbit = leftselected.begin(); mbit != leftselected.end(); ) {
			if (leftsingle.find(mbit->first) != leftsingle.end())
				mbit++;
			else
				leftselected.erase(mbit++);
		}
		
		if (alpha > 1e-8) {
			for (mbit = rightselected.begin(); mbit != rightselected.end(); ) {
				if (mbit->second)
					mbit++;
				else
					rightselected.erase(mbit++);
			}
		}
		else {
			for (mbit = rightselected.begin(); mbit != rightselected.end(); ) {
				if (rightsingle.find(mbit->first) != rightsingle.end())
					mbit++;
				else
					rightselected.erase(mbit++);
			}
		}
	}
}

int BiCover::update(vector<IDPair>& left, vector<IDPair>& right) {
	list<LINT>::iterator lit;
	map<LINT,bool>::iterator mit;
	set<LINT>::iterator iter;
	bool checkside = false;
	
	// init left right labels
	vector<IDPair>().swap(left);
	vector<IDPair>().swap(right);

	// collect left and right labels
	/*
	for (mit = leftselected.begin(); mit != leftselected.end(); mit++) {
		if (mit->second)
			left.push_back(IDPair(mit->first/gsize,mit->first%gsize));
	}
	for (mit = rightselected.begin(); mit != rightselected.end(); mit++) {
		if (mit->second) 
			right.push_back(IDPair(mit->first%gsize,mit->first/gsize));
	}
	*/
	for (mit = leftselected.begin(); mit != leftselected.end(); mit++) {
		if (mit->second) {
			checkside = false;
			for (lit = leftlist[mit->first].begin(); lit != leftlist[mit->first].end(); lit++) {
				if (rightselected[*lit]) {
					checkside = true;
					iter = rightsingle.find(*lit);
					if (iter != rightsingle.end()) {
						right.push_back(IDPair(*lit%gsize,*lit/gsize));
						rightsingle.erase(iter);
					}
				}
			}
			if (checkside)
				left.push_back(IDPair(mit->first/gsize,mit->first%gsize));
			else
				leftsingle.insert(mit->first);
		}
	}
	for (mit = rightselected.begin(); mit != rightselected.end(); mit++) {
		if (mit->second) {
			checkside = false;
 			for (lit = rightlist[mit->first].begin(); lit != rightlist[mit->first].end(); lit++) {
				if (leftselected[*lit]) {
					checkside = true;
					iter = leftsingle.find(*lit);
					if (iter != leftsingle.end()) {
						left.push_back(IDPair(*lit/gsize,*lit%gsize));
						leftsingle.erase(iter);
					}
				}
			} 
			if (checkside)
				right.push_back(IDPair(mit->first%gsize,mit->first/gsize));
			else
				rightsingle.insert(mit->first);
		}
	}

	// update all data structure
	update();

	// for test to check the updated results
	#ifdef DEBUG	
	cout << "------------------------left Labels---------------------------" << endl;
	vector<IDPair>::const_iterator vcit;
	for (vcit = left.begin(); vcit != left.end(); vcit++) {
		cout << "[" << vcit->first << "," << vcit->second << "] ";
	}
	cout << endl;
	cout << "------------------------right Labels---------------------------" << endl;
	for (vcit = right.begin(); vcit != right.end(); vcit++) {
		cout << "[" << vcit->first << "," << vcit->second << "] ";
	}
	cout << endl;
	printBipartiteGraph();
	cout << "------------------------Covered after updated---------------------------" << endl;
	Util::printPairVec(covered); 
	#endif
	
	return leftlist.size()+rightlist.size();
}

// greedily select the object which maximally cover the uncovered pairs
// sidno: 0 left, 1 right, 2 both
pair<LINT,bool> BiCover::selectObj(int sideno) {
	bool leftmark = false;
	LINT sid = MAX_VAL;
	int maxnum = MIN_VAL;
	if (sideno == 0 || sideno == 2) {
		if (leftps.size()>0) {
			leftmark = true;
			sid = leftps.rbegin()->second;
			maxnum = leftps.rbegin()->first;
		}
		else {
			if (sideno != 2)
				cout << "Error selection1" << endl;
		}
	}
	if (sideno == 1 || sideno == 2) {
		if (rightps.size()>0) {
			if (rightps.rbegin()->first>maxnum) {
				leftmark = false;
				sid = rightps.rbegin()->second;
			}
		}
		else {
			if (sideno != 2)
				cout << "Error selection2" << endl;
		}
	}
	
	return pair<LINT,bool>(sid,leftmark);
}

// update the status when object(selectid) is selected
void BiCover::updateState(LINT selectid, bool leftmark) {
	list<LINT>::iterator lit;
	vector<LabelPair>::iterator vit;
	int u,v,deg,leftval,rightval,sc,tc;
	bool leftse, rightse;
	IntPair newpair;
	if (leftmark) {
		#ifdef STAT_ALPHA
		sc = 0;
		tc = 0;
		#endif
		leftselected[selectid] = true;
		for (lit = leftlist[selectid].begin(); lit != leftlist[selectid].end(); lit++) {
			if (!covered[selectid/gsize][(*lit)/gsize]) {
				if (rightselected[*lit]) {
					u = selectid/gsize;
					v = (*lit)/gsize;
					covered[u][v] = true;
					num_covered++;
					
					#ifdef STAT_ALPHA
					rtv[*lit]--;
					#endif
					
					// update leftcover and rightcover, leftps and rightps
					for (vit = spmap[u][v].begin(); vit != spmap[u][v].end(); vit++) {
						leftse = leftselected[vit->first];
						rightse = rightselected[vit->second];
						leftval = 0;
						rightval = 0;
						// identify update value for both objects
						if (!leftse) {
							if (!rightse) { leftval = rightval = (int)(alpha*MFACTOR); }
							else leftval = MFACTOR; 
						}
						else if (!rightse)
							rightval = MFACTOR;
						
						if (leftval>0) {
							deg = leftcover[vit->first];
							leftps.erase(IntPair(deg,(vit->first)));
							if (deg==leftval) {
								leftcover.erase(vit->first);
							}
							else {
								leftcover[vit->first] -= leftval;
								newpair = IntPair(deg-leftval,(vit->first));
								leftps.insert(newpair);
							}
						}
						if (rightval>0) {
							deg = rightcover[vit->second];
							rightps.erase(IntPair(deg,(vit->second)));
							if (deg==rightval) {
								rightcover.erase(vit->second);
							}
							else {
								rightcover[vit->second] -= rightval;
								newpair = IntPair(deg-rightval,(vit->second));
								rightps.insert(newpair);
							}							
						}
					}
				//	spmap[u][v].clear();
					vector<LabelPair>().swap(spmap[u][v]);
					spmap[u].erase(v);
				}
				else {
					// move unselected one to selected part
					if (alpha<1.0) {
						deg = rightcover[*lit];
						rightps.erase(IntPair(deg,(*lit)));
						newpair = IntPair(deg+(int)((1-alpha)*MFACTOR),(*lit));
						rightps.insert(newpair);
						rightcover[*lit] += (int)((1-alpha)*MFACTOR);
					}
					#ifdef STAT_ALPHA
					tc++;
					#endif
				}
				#ifdef STAT_ALPHA
				sc++;
				#endif
			}
			else if (!rightselected[*lit]) {
					// move unselected one to selected part
					if (alpha<1.0) {
						deg = rightcover[*lit];
						rightps.erase(LabelPair(deg,(*lit)));
						newpair = IntPair(deg+(int)((1-alpha)*MFACTOR),(*lit));
						rightps.insert(newpair);
						rightcover[*lit] += (int)((1-alpha)*MFACTOR);
					}
			}
		}
		#ifdef STAT_ALPHA
		lsv[selectid] = sc;
		ltv[selectid] = tc;
		if (tc > sc) cout << "error tc = " << tc << "\tsc = " << sc << endl;
		if (sc > 0) {
			real_alpha2 += (tc*1.0)/(sc*1.0);
			tlsize++;
		}
		#endif
		leftps.erase(IntPair(leftcover[selectid],selectid));
		leftcover.erase(selectid);
	}
	else {
		#ifdef STAT_ALPHA
		sc = 0;
		tc = 0;
		#endif
		rightselected[selectid] = true;
		for (lit = rightlist[selectid].begin(); lit != rightlist[selectid].end(); lit++) {
			if (!covered[*lit/gsize][selectid/gsize]) {
				if (leftselected[*lit]) {
					u = *lit/gsize;
					v = selectid/gsize;
					covered[u][v] = true;
					num_covered++;
					
					#ifdef STAT_ALPHA
					ltv[*lit]--;
					#endif

					// update leftcover and rightcover
					for (vit = spmap[u][v].begin(); vit != spmap[u][v].end(); vit++) {
						leftse = leftselected[vit->first];
						rightse = rightselected[vit->second];
						leftval = 0;
						rightval = 0;
						// identify update value for both objects
						if (!leftse) {
							if (!rightse) { leftval = rightval = (int)(alpha*MFACTOR); }
							else leftval = MFACTOR; 
						}
						else if (!rightse)
							rightval = MFACTOR;
							
						// update leftcover and rightcover
						if (leftval>0) {
							deg = leftcover[vit->first];
							leftps.erase(IntPair(deg,(vit->first)));
							if (deg==leftval) {
								leftcover.erase(vit->first);
							}
							else {
								leftcover[vit->first] -= leftval;
								newpair = IntPair(deg-leftval,(vit->first));
								leftps.insert(newpair);
							}
						}
						if (rightval>0) {
							deg = rightcover[vit->second];
							rightps.erase(IntPair(deg,(vit->second)));
							if (deg==rightval) {
								rightcover.erase(vit->second);
							}
							else {
								rightcover[vit->second] -= rightval;
								newpair = IntPair(deg-rightval,(vit->second));
								rightps.insert(newpair);
							}
						}						
					}
				//	spmap[u][v].clear();
					vector<LabelPair>().swap(spmap[u][v]);
					spmap[u].erase(v);
				}
				else {
					// perform move operation
					if (alpha<1.0) {
						deg = leftcover[*lit];
						leftps.erase(IntPair(deg,(*lit)));
						newpair = IntPair(deg+(int)((1-alpha)*MFACTOR),(*lit));
						leftps.insert(newpair);
						leftcover[*lit] += (int)((1-alpha)*MFACTOR);
					}
					#ifdef STAT_ALPHA
					tc++;
					#endif
				}
				#ifdef STAT_ALPHA
				sc++;
				#endif
			}
			else if (!leftselected[*lit]) {
				// perform move operation
				if (alpha<1.0) {
					deg = leftcover[*lit];
					leftps.erase(IntPair(deg,(*lit)));
					newpair = IntPair(deg+(int)((1-alpha)*MFACTOR),(*lit));
					leftps.insert(newpair);
					leftcover[*lit] += (int)((1-alpha)*MFACTOR);
				}			
			}
		}
		#ifdef STAT_ALPHA
		rsv[selectid] = sc;
		rtv[selectid] = tc;
		if (sc > 0) {
			real_alpha2 += (tc*1.0)/(sc*1.0);
			tlsize++;
		}
		#endif
		rightps.erase(LabelPair(rightcover[selectid],selectid));
		rightcover.erase(selectid);
	}
}	

void BiCover::scp(int cycleid, int end_ind) {
	pair<LINT,bool> obj;
	int iter = 0;

	// compute the sum of uncovered shortest paths
	int num_sp = 0;
	for (int i = 0; i < covered.size(); i++)
		num_sp += covered[i].size();
	double cover_ratio = 1.0;
	/*
	if (cycleid==0 && end_ind<gsize)
		cover_ratio = 0.9;
	*/
		
	num_covered = 0;
	while (num_covered < num_sp*cover_ratio) {
		obj = selectObj(2);

		#ifdef DEBUG
		cout << "Iteration " << iter << " select  ID " << obj.first << " " << obj.second << " (";
		if (obj.second) 
			cout << obj.first/gsize << "," << obj.first%gsize << ") left" << endl;
		else
			cout << obj.first%gsize << "," << obj.first/gsize << ") right" << endl;
		cout << "cost = ";
		if (obj.second)
			cout << leftcover[obj.first];
		else
			cout << rightcover[obj.first];
		cout << "\tnum_sp=" << num_sp;
		cout << "\tnum_covered=" << num_covered << "\tleftps size=" << leftps.size() 
				<< "\trightps size=" << rightps.size() << endl;
		cout << "----------------------------------before update---------------------------" << endl;
		printPairSet(leftps);
		printPairSet(rightps);
		cout << "************************leftcover******************"<< endl;
		Util::printMap(leftcover);
		cout << "************************rightcover*****************"<< endl;
		Util::printMap(rightcover);		
		cout << "------------------------covered------------------------" << endl;
		Util::printPairVec(covered);
		#endif
		
		// update status
		updateState(obj.first,obj.second);
		
		#ifdef DEBUG
		cout << "----------------------------------after update-----------------------------" << endl;
		printPairSet(leftps);
		printPairSet(rightps);
		cout << "************************leftcover******************"<< endl;
		Util::printMap(leftcover);
		cout << "************************rightcover******************"<< endl;
		Util::printMap(rightcover);	
		cout << "------------------------covered------------------------" << endl;
		Util::printPairVec(covered);
		#endif
				
		iter++;
	}
	
	cout << "\t# of iterations in SCP3: " << iter << endl;
}

vector<IDPair> BiCover::leftPairs() {
	list<LINT>::iterator lit;
	map<LINT,bool>::iterator mit;

	vector<IDPair> left;
	// collect left and right labels
	for (mit = leftselected.begin(); mit != leftselected.end(); mit++) {
		if (mit->second) 
			left.push_back(IDPair(mit->first/gsize,mit->first%gsize));
	}
	
	return left;
}
		
vector<IDPair> BiCover::rightPairs() {
	list<LINT>::iterator lit;
	map<LINT,bool>::iterator mit;
	
	vector<IDPair> right;
	for (mit = rightselected.begin(); mit != rightselected.end(); mit++) {
		if (mit->second) 
			right.push_back(IDPair(mit->first%gsize,mit->first/gsize));
	}
	
	return right;
}

double BiCover::stat_alpha() {
	double ra = 0;
	
	#ifdef STAT_ALPHA
	long total_size = 0;
	map<LINT,int>::iterator smit, tmit;
	if (lsv.size() != ltv.size()) {
		cerr << "stat_alpha left side error!" << endl;
		cout << lsv.size() << "\t" << ltv.size() << endl;
		for (smit = lsv.begin(); smit != lsv.end(); smit++)
			cout << smit->first << "|" << smit->second << "\t";
		cout << endl;
		for (smit = ltv.begin(); smit != ltv.end(); smit++)
			cout << smit->first << "|" << smit->second << "\t";
		cout << endl;
	}
	if (rsv.size() != rtv.size()) cerr << "stat_alpha right side error!" << endl;
	for (smit=lsv.begin(), tmit=ltv.begin(); smit != lsv.end(); smit++,tmit++) {
		if (smit->second>0) {
			ra += (tmit->second*1.0)/(smit->second*1.0);
			total_size++;
		}
	}
	for (smit=rsv.begin(), tmit=rtv.begin(); smit != rsv.end(); smit++,tmit++) {
		if (smit->second>0) {
			ra += (tmit->second*1.0)/(smit->second*1.0);
			total_size++;
		}
	}
	ra = ra/(total_size*1.0);
	#endif
	
	return ra;
}

double BiCover::stat_alpha2() {
	int total_size = 0;
	/*
	map<LINT,bool>::const_iterator mcit;
	for (mcit = leftselected.begin(); mcit != leftselected.end(); mcit++)
		if (mcit->second)
			total_size++;
	for (mcit = rightselected.begin(); mcit != rightselected.end(); mcit++)
		if (mcit->second)
			total_size++;		
	*/
	return real_alpha2/(tlsize*1.0);
}

void BiCover::printBipartiteGraph() {
	cout << "*******************BiCover*******************" << endl;
	map<LINT,list<LINT> >::iterator mit;
	list<LINT>::iterator lit;
	for (mit = leftlist.begin(); mit != leftlist.end(); mit++) {
		cout << "(" << mit->first/gsize << "," << mit->first%gsize << ") : ";
		for (lit = mit->second.begin(); lit != mit->second.end(); lit++) {
			cout << "(" << (*lit)%gsize << "," << (*lit)/gsize << ") ";
		}
		cout << endl;
	}
	cout << endl;
	cout << "From right side" << endl;
	for (mit = rightlist.begin(); mit != rightlist.end(); mit++) {
		cout << "(" << mit->first%gsize << "," << mit->first/gsize << ") : ";
		for (lit = mit->second.begin(); lit != mit->second.end(); lit++) {
			cout << "(" << (*lit)/gsize << "," << (*lit)%gsize << ") ";
		}
		cout << endl;
	}
	cout << endl;	
}

void BiCover::printPairSet(const PairSet& ps) {
	cout << "=========================PairSet=======================" << endl;
	PairSet::const_iterator pit;
	for (pit = ps.begin(); pit != ps.end(); pit++) {
		cout << "[v" << pit->second << "," << pit->first << "] ";
	}
	cout << endl;
}
