//this code is to implement the algorithm 4 in Kewang's draft
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct Pair{
	long u;
	long v;
};

struct y_element{
	long u;
	long d;
};
typedef y_element x_element;
struct y_invert_list{
	vector<y_element> l;
	map<long,long> m;
};

struct entry{
	long u;
	long d;
	void *pointer;
};

struct invert{
	vector<entry> l;
};

map<long,invert> global_invert;
map<long,map<long,long>> global_pair;

typedef y_invert_list x_invert_list;

struct in_element{
	long x;
	long d;
};

struct out_elemnt{
	long y;
	long d;
};


struct answer_entry{
	vector<long> e;
};

//node u_k is an id, w_k is a word
vector<long> getInverted(string key)// this function is equal to Pi in the draft
//input is the keyword, output is the list of nodes that contain this keyword
{
	
}

vector<in_element> IN(long u, long d)
{// this function is equal to IN_d(u_k) in kewang's draft

}

vector<out_element> OUT(long u, long d)
{

}

map<long,y_element> constructY(string key,long D)
{
	map<long,y_element> Y;
	vector<long> key_invert_list = getInverted(key);
	for(long i=0; i<key_invert_list.size(); i++)
	{
		long u = key_invert_list[i];
		vector<in_element> in_u = IN(u,D);
		for(long j=0; j<in_u.size(); j++)
		{
			long x = in_u[j].x;
			long d = in_u[j].d;
			if(Y.find(x) == Y.end())
			{
				y_invert_list y;
				y.m[u]=d;
				Y[x] = y;
			}
			else
			{
				map<long,long> m = Y[x].m;
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

long nearest_ancestor(long y)
{
	long x;
	//use the T-structure to find the ancester
	return x;// x is ancester of y;
}

long d_tree(long x, long y)
{
	return 2;
}

bool check(vector<long> a, long v, long D)
{
	for(long i=0; i<a.size(); i++)
	{
		long u = a[i];
		map<long,map<long,long>> ::iterator mmit;
		mmit = global_pair.find(u);
		map<long,long> :: iterator mit;
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

void step2(string key, long D, map<long,long> nodeAnswer, vector<answer_entry> answer )
{
	map<long,y_invert_list> Y = constructY( key, D);
	map<long,x_invert_list> X;
	for(map<long,y_invert_list> it=Y.begin(); it!=Y.end(); it++)
	{
		for(long i=0; i<(*it).second.l.size(); i++)
		{
			long y = (*it).first;
			long u = (*it).second.l[i].u;
			long d = (*it).second.l[i].d;
			long x = nearest_ancestor(y,nodeAnswer);
			long dT =  d_tree(x,y);
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
					map<long,long> m = X[x].m;
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
	for(map<long,x_invert_list>:: iterator it = X.begin(); it!=X.end(); it++)
	{
		long x = (*it).first;
		invert x_x_l = global_invert[x];
		x_invert_list x_y_l = (*it).second;
		for(long i=0; i<x_x_l.size(); i++)
		{
			long uj = x_x_l[i].u;
			long dj = x_x_l[i].d;
			for(long j=0; j<x_y_l.size(); j++)
			{
				long uk = x_y_l[j].u;
				long dk = x_y_l[j].d;
				if( (uj + uk) <D)
				{
					map<long,map<long,long>> ::iterator mmit;
					mmit = global_pair.find(uj);
					map<long,long> :: iterator mit;
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
						map<long,long> mmap;
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
	vector<long> key_list = getInverted( key);
	for(long i=0; i<answer.size(); i++)
	{
		for(long j=0; j<key_list.size(); j++)
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
