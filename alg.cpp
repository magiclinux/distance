//this code is to implement the algorithm 4 in Kewang's draft
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct Pair{
	int u;
	int v;
};

struct y_element{
	int u;
	int d;
};
typedef y_element x_element;
struct y_invert_list{
	vector<y_element> l;
	map<int,int> m;
};

struct entry{
	int u;
	int d;
	void *pointer;
};

struct invert{
	vector<entry> l;
};

map<int,invert> global_invert;
map<int,map<int,int>> global_pair;

typedef y_invert_list x_invert_list;

struct in_element{
	int x;
	int d;
};

struct out_elemnt{
	int y;
	int d;
};


struct answer_entry{
	vector<int> e;
};

//node u_k is an id, w_k is a word
vector<int> getInverted(string key)// this function is equal to Pi in the draft
//input is the keyword, output is the list of nodes that contain this keyword
{
	
}

vector<in_element> IN(int u, int d)
{// this function is equal to IN_d(u_k) in kewang's draft

}

vector<out_element> OUT(int u, int d)
{

}

map<int,y_element> constructY(string key,int D)
{
	map<int,y_element> Y;
	vector<int> key_invert_list = getInverted(key);
	for(int i=0; i<key_invert_list.size(); i++)
	{
		int u = key_invert_list[i];
		vector<in_element> in_u = IN(u,D);
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
	return Y;
}

int nearest_ancestor(int y)
{
	int x;
	//use the T-structure to find the ancester
	return x;// x is ancester of y;
}

int d_tree(int x, int y)
{
	return 2;
}

bool check(vector<int> a, int v, int D)
{
	for(int i=0; i<a.size(); i++)
	{
		int u = a[i];
		map<int,map<int,int>> ::iterator mmit;
		mmit = global_pair.find(u);
		map<int,int> :: iterator mit;
		if(mmit != global_pair.end())
		{
			mit = (*mmit).find(v);
			if(mit != (*mmit).end())
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

void step2(string key, int D, map<int,int> nodeAnswer, vector<answer_entry> answer )
{
	map<int,y_invert_list> Y = constructY( key, D);
	map<int,x_invert_list> X;
	for(map<int,y_invert_list> it=Y.begin(); it!=Y.end(); it++)
	{
		for(int i=0; i<(*it).second.l.size(); i++)
		{
			int y = (*it).first;
			int u = (*it).second.l[i].u;
			int d = (*it).second.l[i].d;
			int x = nearest_ancestor(y,nodeAnswer);
			int dT =  d_tree(x,y);
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
	// use a map to store all the d()
	for(map<int,x_invert_list>:: iterator it = X.begin(); it!=X.end(); it++)
	{
		int x = (*it).first;
		invert x_x_l = global_invert[x];
		x_invert_list x_y_l = (*it).second;
		for(int i=0; i<x_x_l.size(); i++)
		{
			int uj = x_x_l[i].u;
			int dj = x_x_l[i].d;
			for(int j=0; j<x_y_l.size(); j++)
			{
				int uk = x_y_l[j].u;
				int dk = x_y_l[j].d;
				if( (uj + uk) <D)
				{
					map<int,map<int,int>> ::iterator mmit;
					mmit = global_pair.find(uj);
					map<int,int> :: iterator mit;
					if(mmit != global_pair.end())
					{
						mit = (*mmit).find(uk);
						if(mit != (*mmit).end())
						{
							if( (*mit).second > (uj+uk) )
								(*mmit)[uk] = uj+uk;
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
	vector<answer_entry> new_answer;
	vector<int> key_list = getInverted( key);
	for(int i=0; i<answer.size(); i++)
	{
		for(int j=0; j<key_list.size(); j++)
		{
			if(check(answer[i],key_list[j],D))
			{
				answer_entry a_e;
				a_e = answer[i];
				a_e.pysh_back(key_list[j]);
				new_answer.push_back(a_e);
				
			}
		}
		
			
		
	}
}

int main()
{
	
}
