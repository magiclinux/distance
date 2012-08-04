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
	labelfile = "label.txt";
	sptreefile="sptree.txt";
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
	writeSptree(sptree);
	// compute all pairs shortest paths for sptree
	GraphUtil::computeShortestDistanceOnTree(sptree,treedist,gsize);
	ReducedGraph().swap(sptree);
	
	// for test
	//#ifdef DEBUG
	showLabels();
	Util::printSparseVec(treedist);
	//#endif
	
	
	writeLabels();
}

void Distance::writeAll()
{
	
	//Util::writeSparseVec(treedist);
}

void Distance::readAll()
{
	readSptree();
	readLabels();
	//Util::readSparseVec(treedist);
	
	
}

void Distance::writeSptree(const ReducedGraph& rg)
{
	printf("write sptree begin ...\n");
	ReducedGraph::const_iterator rit;
	vector<int>::const_iterator vit;
	
	FILE * fp = fopen(sptreefile.c_str(),"w");
	fprintf(fp,"%d\n",rg.size());
	for (rit = rg.begin(); rit != rg.end(); rit++) {
		fprintf(fp,"%d %d\n", rit->first,rit->second.size());
		for(vit = rit->second.begin(); vit != rit->second.end(); vit++)
			fprintf(fp,"%d ",*vit);
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	printf("write sptree end ...\n");
}

void Distance::readSptree()
{
	printf("read sptree begin ...\n");
	ReducedGraph rg;
	FILE * fp = fopen(sptreefile.c_str(),"r");
	int n,m,id,a;
	fscanf(fp,"%d",&n);
	for(int i=0; i<n; i++)
	{
		fscanf(fp,"%d %d",&id,&m);
		vector<int> v;
		for(int j=0; j<m; j++)
		{
			fscanf(fp,"%d",&a);
			v.push_back(a);
		}
		rg[id] = v;
		
	}
	fclose(fp);
	printf("read sptree end ...\n");
	printf("recover treedist\n");
	GraphUtil::computeShortestDistanceOnTree(rg,treedist,g.num_vertices());// I hope it works.
	printf("recover linearTree\n");
	reducedGraphToLinear(rg,g.num_vertices());
}

map<int,int> Distance::answerToNodeAnswer(vector<vector<int> > answer)
{
	map<int , int > m;
	for(int i=0; i<answer.size(); i++)
	{
		for(int j=0; j<answer[i].size(); j++)
			m[answer[i][j]]=1;
	}
	return m;
}

void Distance::printAnswer(vector<vector<int> > answer)
{
	cout <<"The answer:"<< endl;
	for(int i=0; i<answer.size(); i++)
	{
		cout << "["<<i<<"] ";
		for(int j=0; j<answer[i].size(); j++)
			cout << answer[i][j]<< " ";
		cout <<endl;
	}
	cout << endl;
}

void Distance::printNodeAnswer(map<int,int> nodeAnswer)
{
	cout <<"Print the unique node in the answer"<<endl;
	for(map<int,int> ::iterator it = nodeAnswer.begin(); it!= nodeAnswer.end(); it++)
		cout << (*it).first<<" ";
	cout << endl;
}

void Distance::printY(map<int,y_invert_list > Y)
{
	cout << "Look at Y:"<<endl;
	cout << "Y.size : "<< Y.size()<<endl;
	
	for(map<int,y_invert_list > :: const_iterator it = Y.begin(); it!=Y.end(); it++)
	{
		cout <<"[" << (*it).first<<"] " << "<"<<(*it).second.m.size()<<","<<(*it).second.l.size() <<"> ";
		for(map<int,int > :: const_iterator mit = (*it).second.m.begin(); mit!= (*it).second.m.end(); mit++)
			cout << "("<<(*mit).first<<","<<(*mit).second<<") ";
		cout << endl;
	}
	cout << endl;
}

void Distance::printX(map<int,x_invert_list > X)
{
	cout << "Look at X:"<<endl;
	cout << "X.size : "<< X.size()<<endl;
	
	for(map<int,x_invert_list > :: const_iterator it = X.begin(); it!=X.end(); it++)
	{
		cout <<"[" << (*it).first<<"] " << "<"<<(*it).second.m.size()<<","<<(*it).second.l.size() <<"> ";
		for(map<int,int > :: const_iterator mit = (*it).second.m.begin(); mit!= (*it).second.m.end(); mit++)
			cout << "("<<(*mit).first<<","<<(*mit).second<<") ";
		cout << endl;
	}
	cout << endl;
}

void Distance::printGlobalPair()
{
	cout<< "Print global Pairs : "<<endl;
	for(map<int,map<int,int> > ::iterator mmit= global_pair.begin(); mmit!=global_pair.end(); mmit++)
	{
		for(map<int,int> :: iterator mit = (*mmit).second.begin(); mit!=(*mmit).second.end(); mit++ )
			cout <<	"("<<(*mmit).first<<","<<(*mit).first<<","<<(*mit).second<<") ";
	}
	cout<<endl;
}

void Distance::constructGlobalInvert()
{
	for(int i=0; i< l_out.size(); i++)
	{
		for(map<int,int> :: iterator mit = l_out[i].begin(); mit!=l_out[i].end(); mit++)
		{
			int u=i;
			int x=(*mit).first;
			int d=(*mit).second;
			map<int,vector<entry> > :: iterator git=global_invert.find(x);
			if(git == global_invert.end())
			{
				entry e; e.u=u; e.d=d;
				vector<entry> v; v.push_back(e);
				global_invert[x] = v;
			}
			else
			{
				entry e; e.u=u; e.d=d;
				(*git).second.push_back(e);
			}
		}
	}
}

void Distance::printGlobalInvert()
{
	cout << "Print the global Invert : "<<endl;
	cout << "size : "<<global_invert.size()<<endl;
	for(map<int,vector<entry> > :: iterator it= global_invert.begin(); it!=global_invert.end(); it++)
	{
		cout << "["<<(*it).first<<"] ";
		for(int i=0; i< (*it).second.size(); i++)
			cout << "("<<(*it).second[i].u << "," << (*it).second[i].d << ") ";
		cout << endl;
	}
}

void Distance::run(vector<string> keys,int D)
{

	map<int,int> nodeAnswer;
	vector<vector<int> > answer;
	vector<int> l = g.getKeyword(keys[0]);
	for(int i=0; i<l.size(); i++)
	{
		vector<int> v;
		v.push_back(l[i]);
		answer.push_back(v);
	}
	constructGlobalInvert();
	//printGlobalInvert();
	for(int i=1; i<keys.size(); i++)
	{
		//printAnswer(answer);
		nodeAnswer =  answerToNodeAnswer(answer);
		//printNodeAnswer(nodeAnswer);
		answer = step2(keys[i],D,nodeAnswer,answer);
		///printAnswer(answer);
		
	}
	printAnswer(answer);
}

vector<in_element> Distance::getIn(int u, int D)
{
	vector<in_element> in_list;
	map<int,int> ::const_iterator it;
	
	for(it= l_in[u].begin(); it!= l_in[u].end(); it++)
	{
		in_element ie;
		ie.x = (*it).first;
		ie.d = (*it).second;
		if(ie.d <= D)
		in_list.push_back(ie);
	}
	return in_list;
}
vector<out_element> Distance::getOut(int u, int D)
{
	vector<out_element> out_list;
	map<int,int> ::const_iterator it;
	
	for(it= l_out[u].begin(); it!= l_out[u].end(); it++)
	{
		out_element oe;
		oe.y = (*it).first;
		oe.d = (*it).second;
		if(oe.d <= D)
		out_list.push_back(oe);
	}
	return out_list;
}

map<int,y_invert_list> Distance::constructY(string key,int D)
{
	map<int,y_invert_list > Y;
	vector<int> key_invert_list = g.getKeyword(key);
	for(int i=0; i<key_invert_list.size(); i++)
	{
		int u = key_invert_list[i];
		vector<in_element > in_u = getIn(u,D);
		for(int j=0; j<in_u.size(); j++)
		{
			int x = in_u[j].x;
			int d = in_u[j].d;
			if(Y.find(x) == Y.end())
			{
				y_invert_list y;
				y.m[u]=d;
				Y[x] = y;
			}
			else
			{
				map<int,int> m = Y[x].m;
				if(m.find(u)==m.end())
				{
					m[u] = d;
					Y[x].m = m;
					
				}
				else if(m[u] > d)
				{
					m[u] = d;
					Y[x].m = m;// we also need to convert all the things from map to list and sort them
				}
				
			}
		}
		
	}


	//push the data from map to vector

	for(map<int,y_invert_list > :: iterator it = Y.begin(); it!=Y.end(); it++)
	{
		
		for(map<int,int > :: iterator mit = (*it).second.m.begin(); mit!= (*it).second.m.end(); mit++)
		{
			y_element ye;
			ye.u = (*mit).first; ye.d =(*mit).second;
			(*it).second.l.push_back(ye);
		}
		
	}
	// print the data Y
	//printY(Y);
	return Y;
}

int Distance::nearestAncestor(int y, map<int,int> refs)
{
	int x=y;
	while(x!=-1 && linearTree[x]!=-1)
	{
		x = linearTree[x];
		if(refs.find(x)!=refs.end())
		{
			return x;// x is ancester of y;
		}
	}
	return -1;// there is no ancestor of y belongs to refs, what should I return? Ask kewang tomorrow.
}


bool Distance::check(vector<int> a, int v, int D)
{
	for(int i=0; i<a.size(); i++)
	{
		int u = a[i];
		if(u == v)
			continue;
		map<int,map<int,int> > ::iterator mmit;
		mmit = global_pair.find(u);
		map<int,int> :: iterator mit;
		if(mmit != global_pair.end())
		{
			mit = (*mmit).second.find(v);
			if(mit != (*mmit).second.end())
			{
				if( (*mit).second > D )
					return false;
			}
			else
				return false;

		}
		else// uj doesn't exist in the map
		{
			return false;
		}
	}
	return true;
}



vector<vector<int> > Distance::step2(string key, int D, map<int,int> nodeAnswer, vector<vector<int> > answer )
{
	map<int,y_invert_list> Y = constructY( key, D);
	map<int,x_invert_list> X;
	for(map<int,y_invert_list>::iterator it=Y.begin(); it!=Y.end(); it++)
	{
		for(int i=0; i<(*it).second.l.size(); i++)
		{
			int y = (*it).first;
			int u = (*it).second.l[i].u;
			int d = (*it).second.l[i].d;
			int x = nearestAncestor(y,nodeAnswer);
			if(x == -1)
			{
				cout <<"Error: x is -1"<<endl;
				continue;
			}
				
			int dT =  distance(x,y);
			if( ( d + dT ) <= D )
			{
				if(X.find(x) == X.end())
				{
					x_invert_list x_i;
					x_i.m[u] = d;
					X[x]= x_i;
				}
				else 
				{
					map<int,int> m = X[x].m;
					if(m.find(u) == m.end())
					{
						m[u] = d;
						X[x].m = m;
					}
					else if(m[u] > ( d + dT ) )
					{
						m[u] = d;
						X[x].m = m;
					}
				}
			}
		}
		
	}

	//push the data from map to vector

	for(map<int,x_invert_list > :: iterator it = X.begin(); it!=X.end(); it++)
	{
		
		for(map<int,int > :: iterator mit = (*it).second.m.begin(); mit!= (*it).second.m.end(); mit++)
		{
			x_element xe;
			xe.u = (*mit).first; xe.d =(*mit).second;
			(*it).second.l.push_back(xe);
		}
		
	}
	//print the data X
	//printX(X);

	//printGloablInvert
	//printGlobalInvert();
	
	// use a map to store all the d()
	for(map<int,x_invert_list>:: iterator it = X.begin(); it!=X.end(); it++)
	{
		int x = (*it).first;
		vector<entry> x_x_l = global_invert[x];
		vector<y_element> x_y_l = (*it).second.l;
		if(x_x_l.size()==0)
			cout << "x_x_l size is zero"<<endl;
		for(int i=0; i<x_x_l.size(); i++)
		{
			int uj = x_x_l[i].u;
			int dj = x_x_l[i].d;
			for(int j=0; j<x_y_l.size(); j++)
			{
				int uk = x_y_l[j].u;
				int dk = x_y_l[j].d;
				cout <<"ui,uk "<<uj<<" "<<dj<< " "<<uk<<" "<<dk<<endl;
				if( (dj + dk) <D)
				{
					map<int,map<int,int> > ::iterator mmit;
					mmit = global_pair.find(uj);
					map<int,int> :: iterator mit;
					if(mmit != global_pair.end())
					{
						mit = (*mmit).second.find(uk);
						if(mit != (*mmit).second.end())
						{
							if( (*mit).second > (dj+dk) )
								(*mmit).second[uk] = dj+dk;
						}
					}
					else// uj doesn't exist in the map
					{
						map<int,int> mmap;
						mmap[uk] = dj + dk;
						global_pair[uj] = mmap;
					}
					
				}
				else
				{
					break;
				}
			}
		}
	}

	//print golbal pairs
	//printGlobalPair();

	
	vector<vector<int> > new_answer;
	vector<int> key_list = g.getKeyword( key);
	for(int i=0; i<answer.size(); i++)
	{
		for(int j=0; j<key_list.size(); j++)
		{
			if(check(answer[i],key_list[j],D))
			{
				vector<int> a_e;
				a_e = answer[i];
				a_e.push_back(key_list[j]);
				new_answer.push_back(a_e);
				
			}
		}
		
			
		
	}
	return new_answer;
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

void Distance::writeLabels()
{
	printf("write the labels begin ...\n");
	FILE * fp = fopen(labelfile.c_str(),"w");
	fprintf(fp,"%d\n",l_out.size());
	for(int i=0; i<l_out.size(); i++)
	{
		fprintf(fp,"%d %d\n",i,l_out[i].size());
		for(map<int,int>::iterator mit = l_out[i].begin(); mit != l_out[i].end(); mit++)
			fprintf(fp,"%d %d ",mit->first,mit->second);
		fprintf(fp,"\n");
	}
	
	fprintf(fp,"%d\n",l_in.size());
	for(int i=0; i<l_in.size(); i++)
	{
		fprintf(fp,"%d %d\n",i,l_in[i].size());
		for(map<int,int>::iterator mit = l_in[i].begin(); mit != l_in[i].end(); mit++)
			fprintf(fp,"%d %d ",mit->first,mit->second);
		fprintf(fp,"\n");
	}
	fclose(fp);
	printf("write the labels end ...\n");
}

void Distance::readLabels()
{
	printf("read the labels begin ...\n");
	FILE * fp = fopen(labelfile.c_str(),"r");
	int n;
	fscanf(fp,"%d",&n);
	for(int i=0; i<n; i++)
	{
		int id,m;
		fscanf(fp,"%d%d",&id,&m);
		map<int,int> amap;
		int a,b;
		for(int j=0; j<m; j++)
		{
			fscanf(fp,"%d%d",&a,&b);
			amap[a]=b;
		}
		l_out.push_back(amap);
	}
	
	fscanf(fp,"%d",&n);
	for(int i=0; i<n; i++)
	{
		int id,m;
		fscanf(fp,"%d%d",&id,&m);
		map<int,int> amap;
		int a,b;
		for(int j=0; j<m; j++)
		{
			fscanf(fp,"%d%d",&a,&b);
			amap[a]=b;
		}
		l_in.push_back(amap);
	}
	printf("write the labels end ...\n");
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

