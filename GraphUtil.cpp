#include "GraphUtil.h"

// build maximal weighted spanning tree and rank the vertices based on tc size of their BFS tree
void GraphUtil::buildMaxWST(Graph& g, ReducedGraph& sptree, vector<int>& rank) {
	ReducedGraph::iterator rit;
	vector<int>::iterator vit;
	map<int,int>::iterator mit;
	SparseVec::iterator svit;
	int gsize = g.num_vertices();
	map<int,pair<int,int> > pred_succ;
	SparseVec weight = SparseVec(gsize,map<int,int>());
	multimap<int,int> mrank;
	int tree_tc = 0;

	for (int i = 0; i < gsize; i++) {
		ReducedGraph tree;
		tree_tc = BFSTree(g, i, tree);
		#ifdef ALLPATHS
		treestat(tree, pred_succ, i, gsize);
		#endif
		for (rit = tree.begin(); rit != tree.end(); rit++) {
			for (vit = rit->second.begin(); vit != rit->second.end(); vit++)
				#ifdef ALLPATHS
				weight[rit->first][*vit] += pred_succ[rit->first].first*pred_succ[*vit].second;
				#else
				weight[rit->first][*vit] += 1;
				#endif
		}
		#ifdef ALLPATHS
		map<int,pair<int,int> >().swap(pred_succ);
		#endif
		ReducedGraph().swap(tree);
		mrank.insert(make_pair(tree_tc,i));
	}
	
	// set ranking results
//	rank.clear();
	vector<int>().swap(rank);
	multimap<int,int>::const_reverse_iterator mcit;
	for (mcit = mrank.rbegin(); mcit != mrank.rend(); mcit++)
		rank.push_back(mcit->second);
//	mrank.clear();
	multimap<int,int>().swap(mrank);

	multimap<int,pair<int,int> > wmap;
	multimap<int,pair<int,int> >::reverse_iterator mmit; 
	for (int i = 0; i < weight.size(); i++) {
		for (mit = weight[i].begin(); mit != weight[i].end(); mit++) {
			wmap.insert(make_pair(mit->second,make_pair(i,mit->first)));
		}
	}
//	weight.clear();
	SparseVec().swap(weight);

	// greedy algorithm to extract maximal weight spanning tree
	int numedge = 0;
	Graph tree = Graph(gsize);
	pair<int,int> edge;
	for (mmit = wmap.rbegin(); mmit != wmap.rend(); mmit++) {
		edge = mmit->second;
		if (tree.in_edges(edge.second).size()>0 || checkPath(tree,edge.second,edge.first))
			continue;
		tree.addEdge(edge.first,edge.second);
		numedge++;
		if (numedge>=(gsize-1)) break;
	}

	// compute all pairs shortest paths in the spanning tree
//	computeShortestDistance(tree,treedist);
	
	EdgeList el;
	for (int i = 0; i < gsize; i++) {
		el = tree.out_edges(i);
		if (el.size()>0)
			sptree[i] = vector<int>(el.begin(),el.end());
//		sptree.push_back(vector<int>(el.begin(),el.end()));
//		sptree.insert(el.begin(),el.end());
	}
}

void GraphUtil::buildSPTree(Graph& g, ReducedGraph& sptree, vector<int>& rank) {
	int max_cn = 0, id, cn, root;
	ReducedGraph::iterator rit;
	hash_map<int,int> order, min;
	multimap<int,int> mrank;
	map<int,pair<int,int> > pred_succ;
	int tree_tc = 0;
	
	int gsize = g.num_vertices();
	for (int i = 0; i < gsize; i++) {
		ReducedGraph tree;
		tree_tc = BFSTree(g, i, tree);
		treestat(tree, pred_succ, i, gsize);
		cn = 0;
		for (rit = tree.begin(); rit != tree.end(); rit++) {
			cn += pred_succ[rit->first].second;
		}
		if (max_cn < cn) {
			sptree.clear();
			sptree = tree;
			max_cn = cn;
		}
		mrank.insert(make_pair(cn,i));
	}
	vector<int>().swap(rank);
	multimap<int,int>::const_reverse_iterator mcit;
	for (mcit = mrank.rbegin(); mcit != mrank.rend(); mcit++)
		rank.push_back(mcit->second);
	multimap<int,int>().swap(mrank);
}

int GraphUtil::genCandLabels(Graph& g, int vid, ReducedGraph& sptree, 
		map<int,set<int> >& tpairs, vector<int>& left, map<int,set<int> >& right) {
	int u;
	deque<int> que;
	vector<int>::const_iterator sit;
	EdgeList el;
	EdgeList::iterator eit;
	ReducedGraph temp_right;
	ReducedGraph treeedges;
	int gsize = g.num_vertices();
	vector<bool> in = vector<bool>(gsize,false);
	vector<bool> visit = vector<bool>(gsize,false);

//	left.clear();
	vector<int>().swap(left);
	left.push_back(vid);
	que.push_back(vid);
	in[vid] = true;
	while (!que.empty()) {
		u = que.front();
		que.pop_front();
		if (!visit[u]) {
			visit[u] = true;
			// visit tree edges firstly
			if (sptree.find(u) != sptree.end()) {
				for (sit = sptree[u].begin(); sit != sptree[u].end(); sit++) {
					if (!visit[*sit] && !in[*sit]) {
						que.push_back(*sit);
						in[*sit] = true;
						treeedges[u].push_back(*sit);
						temp_right[*sit] = temp_right[u];
					}
				}
			}
			// visit regular edges in original graph
			el = g.out_edges(u);
			for (eit = el.begin(); eit != el.end(); eit++) {
				if (!visit[*eit] && !in[*eit]) {
					que.push_back(*eit);
					in[*eit] = true;
					left.push_back(*eit);
					temp_right[*eit] = temp_right[u];
					temp_right[*eit].push_back(u);
				}
			}
		}
	}	
	
	// for test
/* 	
	cout << "--------------------tree edges --------------------" << endl;
	Util::printReducedGraph(treeedges); */
	
	// compute tree TC
//	tpairs.clear();
	map<int,set<int> >().swap(tpairs);
	treetc(treeedges, tpairs);
//	treeedges.clear();
	ReducedGraph().swap(treeedges);

	// include themselves
	ReducedGraph::iterator rit;
	vector<int>::iterator vit;
	for (rit = temp_right.begin(); rit != temp_right.end(); rit++) {
		for (vit = rit->second.begin(); vit != rit->second.end(); vit++)
			right[*vit].insert(rit->first);
		right[rit->first].insert(rit->first);
	}
	
	// stat size
	int size = 0;
	map<int,set<int> >::const_iterator msit;
	for (msit = right.begin(); msit != right.end(); msit++)
		size += msit->second.size();
	size += left.size();
	
	return size;
}

void GraphUtil::computePartialSD(Graph& g, vector<int>& rank, int start_ind, 
				int& end_ind, SparseVec& pdist, long tc_limit) {
	int u;
	deque<int> que;
	EdgeList el;
	EdgeList::iterator eit;
	map<int,int>::iterator mit;

	pdist.clear();
	int gsize = g.num_vertices();
	vector<bool> in = vector<bool>(gsize,false);
	vector<bool> visit = vector<bool>(gsize,false);
	int inc_count = 0;
	int sid = start_ind;
	int i = 0, j = rank[sid];
	while (true) {
		pdist.push_back(map<int,int>());
		pdist[i][j] = 0;
		que.clear();
		que.push_back(j);
		in[j] = true;
		while (!que.empty()) {
			u = que.front();
			que.pop_front();
			if (!visit[u]) {
				visit[u] = true;
				el = g.out_edges(u);
				for (eit = el.begin(); eit != el.end(); eit++) {
					if (!visit[*eit] && !in[*eit]) {
						que.push_back(*eit);
						in[*eit] = true;
						pdist[i][*eit] = pdist[i][u]+1;	
					}
				}
			}
		}
		for (mit = pdist[i].begin(); mit != pdist[i].end(); mit++) {
			in[mit->first] = false;
			visit[mit->first] = false;
		}
		i++;
		sid++;
		if (sid < rank.size()) 
			j = rank[sid];
		inc_count += pdist[i].size();
		
		if (inc_count >= tc_limit)
			break;
	}
	end_ind = sid;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GraphUtil::computeShortestDistanceOnTree(ReducedGraph& sptree, SparseVec& dist, int gsize) {
	int u;
	deque<int> que;
	vector<int>::const_iterator eit;
	map<int,int>::iterator mit;

	dist = SparseVec(gsize,map<int,int>());
	vector<bool> in = vector<bool>(gsize,false);
	vector<bool> visit = vector<bool>(gsize,false);
	for (int i = 0; i < gsize; i++) {
		dist[i][i] = 0;
		que.clear();
		que.push_back(i);
		in[i] = true;
		while (!que.empty()) {
			u = que.front();
			que.pop_front();
			if (!visit[u]) {
				visit[u] = true;
				if (sptree.find(u) != sptree.end()) {
					for (eit = sptree[u].begin(); eit != sptree[u].end(); eit++) {
						if (!visit[*eit] && !in[*eit]) {
							que.push_back(*eit);
							in[*eit] = true;
							dist[i][*eit] = dist[i][u]+1;	
						}
					}
				}
			}
		}
		for (mit = dist[i].begin(); mit != dist[i].end(); mit++) {
			in[mit->first] = false;
			visit[mit->first] = false;
		}
	}
}

void GraphUtil::computeShortestDistance(Graph& g, SparseVec& dist) {
	int u;
	deque<int> que;
	EdgeList el;
	EdgeList::iterator eit;
	map<int,int>::iterator mit;

	int gsize = g.num_vertices();
	dist = SparseVec(gsize,map<int,int>());
	vector<bool> in = vector<bool>(gsize,false);
	vector<bool> visit = vector<bool>(gsize,false);
	for (int i = 0; i < gsize; i++) {
		dist[i][i] = 0;
		que.clear();
		que.push_back(i);
		in[i] = true;
		while (!que.empty()) {
			u = que.front();
			que.pop_front();
			if (!visit[u]) {
				visit[u] = true;
				el = g.out_edges(u);
				for (eit = el.begin(); eit != el.end(); eit++) {
					if (!visit[*eit] && !in[*eit]) {
						que.push_back(*eit);
						in[*eit] = true;
						dist[i][*eit] = dist[i][u]+1;	
					}
				}
			}
		}
		for (mit = dist[i].begin(); mit != dist[i].end(); mit++) {
			in[mit->first] = false;
			visit[mit->first] = false;
		}
	}
}

//void GraphUtil::buildSpanningTree1(Graph& g, map<int,set<int> >& tpairs) {
int GraphUtil::buildSpanningTree1(Graph& g, SparseVec& treedist, SparseVec& dist, ReducedGraph& sptree) { 
	int max_cn = 0, id, cn, root;
	ReducedGraph *max_tree = new ReducedGraph();
	ReducedGraph *temp_tree = new ReducedGraph();
	ReducedGraph::iterator rit;
	hash_map<int,int> order, min;
	int gsize = g.num_vertices();
	for (int i = 0; i < gsize; i++) {
		if (g.out_edges(i).size() <= 0) continue;
		id = cn = 0;
		order.clear();
		min.clear();
		ReducedGraph *tree = new ReducedGraph();
		BFSTree(g, i, *tree);
		PVisit(*tree, i, id, cn, order, min);
		// for test
		/*
		Util::printReducedGraph(*tree);
		cout << "cn = " << cn << endl;
		*/
		if (max_cn < cn) {
			root = i;
			max_cn = cn;
			temp_tree = max_tree;
			max_tree = tree;
			delete temp_tree;
		}
	}
	distPairs(*max_tree, root, treedist, dist);
	sptree = *max_tree;
	/*
	sptree = vector<vector<int> >(gsize,vector<int>());
	for (rit = max_tree->begin(); rit != max_tree->end(); rit++)
		sptree[rit->first] = rit->second;
		*/
//	distPairs(max_tree, root, tpairs);
	
	// for test
//	Util::printReducedGraph(*max_tree);

	//	Util::printSparseVec(treedist);
	/*
	cout << "root = " << root << endl;
	map<int,set<int> >::iterator mit;
	set<int>::iterator sit;
	for (mit = tpairs.begin(); mit != tpairs.end(); mit++) {
		cout << "from " << mit->first << ":";
		for (sit = mit->second.begin(); sit != mit->second.end(); sit++)
			cout << *sit << " ";
		cout << endl;
	}
	*/
	
	int tnum = 0;
	for (int i = 0; i < gsize; i++)
		tnum += treedist[i].size();
	return tnum;
}


int GraphUtil::buildSpanningTree2(Graph& g, SparseVec& treedist, SparseVec& dist, ReducedGraph& sptree) {
	ReducedGraph::iterator rit;
	vector<int>::iterator vit;
	map<int,int>::iterator mit;
	SparseVec::iterator svit;
	int gsize = g.num_vertices();
	SparseVec weight = SparseVec(gsize,map<int,int>());

	for (int i = 0; i < gsize; i++) {
		ReducedGraph tree;
		BFSTree(g, i, tree);
		for (rit = tree.begin(); rit != tree.end(); rit++) {
			for (vit = rit->second.begin(); vit != rit->second.end(); vit++)
				weight[rit->first][*vit] += 1;
		}
		tree.clear();
	}	
	// for test
	/*
	cout << "weight summary ";
	Util::printSparseVec(weight);
	*/

	multimap<int,pair<int,int> > wmap;
	multimap<int,pair<int,int> >::reverse_iterator mmit; 
	for (int i = 0; i < weight.size(); i++) {
		for (mit = weight[i].begin(); mit != weight[i].end(); mit++) {
			wmap.insert(make_pair(mit->second,make_pair(i,mit->first)));
		}
	}
	weight.clear();

	// greedy algorithm to extract maximal weight spanning tree
	int numedge = 0;
	Graph tree = Graph(gsize);
	pair<int,int> edge;
	for (mmit = wmap.rbegin(); mmit != wmap.rend(); mmit++) {
		// for test
//		cout << "(" << mmit->second.first << "," << mmit->second.second << ")" << mmit->first << endl;	
		edge = mmit->second;
		if (tree.in_edges(edge.second).size()>0 || checkPath(tree,edge.second,edge.first))
			continue;
		tree.addEdge(edge.first,edge.second);
		numedge++;
		if (numedge>=(gsize-1)) break;
	}
	// for test
//	tree.writeGraph(cout);

	computeShortestDistance(tree,treedist);
	
	int tnum = 0;
	EdgeList el;
	for (int i = 0; i < gsize; i++) {
		treedist[i].erase(i);
		tnum += treedist[i].size();
		for (mit = treedist[i].begin(); mit != treedist[i].end(); ) {
			if (mit->second != dist[i][mit->first]) {
				treedist[i].erase(mit++);
			}
			else
				mit++;
		}
		el = tree.out_edges(i);
		sptree[i] = vector<int>(el.begin(),el.end());
	//	sptree.push_back(vector<int>(el.begin(),el.end()));
	}

	// for test
//	Util::printSparseVec(treedist);

	return tnum;
}	

bool GraphUtil::checkPath(Graph& g, int sid, int tid) {
	if (sid == tid) return true;
	EdgeList& el = g.out_edges(sid);
	EdgeList::iterator eit;
	for (eit = el.begin(); eit != el.end(); eit++) {
		if (checkPath(g,*eit,tid))
			return true;
	}
	return false;
}


int GraphUtil::BFSTree(Graph& g, int v, ReducedGraph& tree) {
	int u;
	EdgeList el;
	EdgeList::iterator eit;
	int gsize = g.num_vertices();
	vector<bool> visit = vector<bool>(gsize,false);
	vector<bool> in = vector<bool>(gsize,false);
	deque<int> que;
	que.push_back(v);
	in[v] = true;
	tree.insert(make_pair(v,vector<int>()));
	while (!que.empty()) {
		u = que.front();
		que.pop_front();
		if (!visit[u]) {
			visit[u] = true;
			el = g.out_edges(u);
			for (eit = el.begin(); eit != el.end(); eit++) {
				if (!visit[*eit] && !in[*eit]) {
					que.push_back(*eit);
					in[*eit] = true;
					tree[u].push_back(*eit);
				}
			}
		}
	}
	
	int tcsize = 0;
	ReducedGraph::const_iterator rit;
	for (rit = tree.begin(); rit != tree.end(); rit++)
		tcsize += rit->second.size();
		
	return tcsize;
}

void GraphUtil::treestat(ReducedGraph& tree, map<int,pair<int,int> >& pred_succ, int root, int gsize) {
	map<int,set<int> > tpairs;
	map<int,set<int> >::iterator msit;
	set<int>::iterator sit;
	SparseVec treedist;
	ReducedGraph::iterator rit;
	vector<int>::iterator vit;
	
	set<int> vertices;
	for (rit = tree.begin(); rit != tree.end(); rit++) {
		vertices.insert(rit->first);
		for (vit = rit->second.begin(); vit != rit->second.end(); vit++)
			vertices.insert(*vit);
	}
//	Util::printReducedGraph(tree);
	computeShortestDistanceOnTree(tree, treedist, gsize);	
//	Util::printSparseVec(treedist);
	treetc(tree, tpairs);
	/*
	for (msit = tpairs.begin(); msit != tpairs.end(); msit++) {
		cout << msit->first << ": ";
		for (sit = msit->second.begin(); sit != msit->second.end(); sit++)
			cout << *sit << " ";
		cout << endl;
	}
	*/
	for (sit = vertices.begin(); sit != vertices.end(); sit++) {
		pred_succ[*sit] = make_pair(treedist[root][*sit]+1,tpairs[*sit].size()+1);
	}

	// for test	
	/*
	map<int,pair<int,int> >::iterator mit;
	for (mit = pred_succ.begin(); mit != pred_succ.end(); mit++) {
		cout << "v " << mit->first << " [" << mit->second.first << "," << mit->second.second << "]" << endl;
	}
	cout << endl;
	*/
//	exit(0);
}

void GraphUtil::treetc(ReducedGraph& tree, map<int,set<int> >& tpairs) {
	map<int,bool> mark;
	ReducedGraph::iterator rit;
	
	for (rit = tree.begin(); rit != tree.end(); rit++)
		mark[rit->first] = false;
	for (rit = tree.begin(); rit != tree.end(); rit++) {
		if (!mark[rit->first])
			treetcvisit(tree,rit->first,tpairs,mark);
	}
}

void GraphUtil::treetcvisit(ReducedGraph& tree, int v, map<int,set<int> >& tpairs, map<int,bool>& mark) {
	vector<int>::iterator vit;
	mark[v] = true;
	for (vit = tree[v].begin(); vit != tree[v].end(); vit++) {
		mark[*vit] = true;
		treetcvisit(tree, *vit, tpairs, mark);
	}
	if (tree[v].size() > 0) { 
		for (vit = tree[v].begin(); vit != tree[v].end(); vit++) {
			if (tree[*vit].size()>0)
				tpairs[v].insert(tpairs[*vit].begin(),tpairs[*vit].end());
			tpairs[v].insert(*vit);
		}
	}
}

void GraphUtil::treetc1(ReducedGraph& tree, int v, map<int,set<int> >& tpairs) {
	vector<int>::iterator vit;
	for (vit = tree[v].begin(); vit != tree[v].end(); vit++) {
		treetc1(tree, *vit, tpairs);
	}
	if (tree[v].size() > 0) { 
		for (vit = tree[v].begin(); vit != tree[v].end(); vit++) {
			if (tree[*vit].size()>0)
				tpairs[v].insert(tpairs[*vit].begin(),tpairs[*vit].end());
			tpairs[v].insert(*vit);
		}
	}
}

void GraphUtil::distPairs(ReducedGraph& tree, int v, SparseVec& treedist, SparseVec& dist) { 
	map<int,set<int> > tpairs;
	map<int,set<int> >::iterator mit;
	set<int>::iterator sit;
	
	treetc1(tree, v, tpairs);
	treedist = SparseVec(dist.size(),map<int,int>());
	for (mit = tpairs.begin(); mit != tpairs.end(); mit++) {
		for (sit = mit->second.begin(); sit != mit->second.end(); sit++) {
			treedist[mit->first][*sit] = dist[mit->first][*sit];
		}
	}
}

// return the number of shortest paths covered by tree (cn)
void GraphUtil::PVisit(ReducedGraph& tree, int v, int& id, int& cn, hash_map<int,int>& order, hash_map<int,int>& min) {
	vector<int>::iterator vit;
	for (vit = tree[v].begin(); vit != tree[v].end(); vit++) {
		PVisit(tree, *vit, id, cn, order, min);
		order[*vit] = id;
		id++;
	}
	order[v] = id;
	if (tree[v].size() == 0) 
		min[v] = order[v];
	else {
		min[v] = MAX_VAL;
		for (vit = tree[v].begin(); vit != tree[v].end(); vit++) {
			if (min[*vit]<min[v])
				min[v] = min[*vit];
		}
		cn += order[v]-min[v];
	}
}

int GraphUtil::BFSDist(Graph& g, int sid, int tid) {
	if (sid == tid) return 0;
	
	int u;
	EdgeList el;
	EdgeList::iterator eit;
	int gsize = g.num_vertices();
	vector<bool> visit = vector<bool>(gsize,false);
	vector<bool> in = vector<bool>(gsize,false);
	vector<int> sdist = vector<int>(gsize,0);
	deque<int> que;
	que.push_back(sid);
	in[sid] = true;
	sdist[sid] = 0;
	while (!que.empty()) {
		u = que.front();
		que.pop_front();
		if (!visit[u]) {
			visit[u] = true;
			el = g.out_edges(u);
			for (eit = el.begin(); eit != el.end(); eit++) {
				if (!visit[*eit] && !in[*eit]) {
					sdist[*eit] = sdist[u]+1;
					if (*eit == tid) return sdist[*eit];	
					que.push_back(*eit);
					in[*eit] = true;
				}
			}
		}
	}
	return -1; //unreachable
}

// for test 
void GraphUtil::genRandomRatioDAG(int n, double c, char* filename, double ratio) {
	int threshold = (int)(c*200*ratio);
	int cth = (int)(c*200*(1-ratio));
	Graph g;
	int i, j;
	int rand_num;
	for (i = 0; i < n; i++) 
		g.addVertex(i);

	srand(time(NULL));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			rand_num = rand()%(100*n);
			if (j < i) {
				if (rand_num < threshold)
					g.addEdge(i,j);
			}
			else if (j > i) {
				if (rand_num < cth)
					g.addEdge(i,j);
			}
		}
	}

	ofstream out(filename);
	g.writeGraph(out);
}

// for test
void GraphUtil::genRandomDAG(int n, double c, char* filename) {
	int threshold = (int)(c*200);
	Graph g;
	int i, j;
	int rand_num;
	for (i = 0; i < n; i++) 
		g.addVertex(i);

	srand(time(NULL));
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++)
			if (i != j) {
				rand_num = rand()%(100*n);
				if (rand_num < threshold)
					g.addEdge(i,j);
			}
	}

	ofstream out(filename);
	g.writeGraph(out);
}

void GraphUtil::genRandomGraph(int n, double c, char* filename) {
	int threshold = (int)(c*10);
	Graph g;
	int i, j;
	int rand_num;
	for (i = 0; i < n; i++) 
		g.addVertex(i);

	srand(time(NULL));
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			if (i != j) {
				rand_num = rand()%(10*n);
				if (rand_num < threshold)
					g.addEdge(i,j);
			}
	}

	ofstream out(filename);
	g.writeGraph(out);
}

