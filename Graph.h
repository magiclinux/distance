#ifndef _GRAPH_H
#define _GRAPH_H

#include "Util.h"

using namespace std;
using namespace __gnu_cxx;

#define MAX_VAL 10000000
#define MIN_VAL -10000000

/* struct Vertex {
	bool visited;
}; */

typedef int Vertex;
typedef vector<int> EdgeList;	// edge list represented by vertex id list
typedef vector<Vertex> VertexList;	// vertices list (store real vertex property) indexing by id

struct In_OutList {
	EdgeList inList;
	EdgeList outList;
};
typedef vector<In_OutList> GRA;	// index graph

class Graph {
	private:
		GRA graph;
		VertexList vl;
		int vsize;
		map< string, vector<int> > keyword;
	public:
		Graph();
		Graph(int);
		Graph(GRA&, VertexList&);
		Graph(istream&);
		~Graph();
		void readGraph(istream&);
		void readKeyword(istream&);
		vector<int> getKeyword(string key);
		void writeGraph(ostream&);
		void printGraph();
		void addVertex(int);
		void addEdge(int, int);
		int num_vertices();
		int num_edges();
		VertexList& vertices();
		EdgeList& out_edges(int);
		EdgeList& in_edges(int);
		int out_degree(int);
		int in_degree(int);
		vector<int> getRoots();
		bool hasEdge(int, int);	
		Graph& operator=(const Graph&);
		Vertex& operator[](const int);
		
		void clear();
		void strTrimRight(string& str);

		Graph(hash_map<int,vector<int> >& inlist, hash_map<int,vector<int> >& outlist);
		void extract(hash_map<int,vector<int> >& inlist, hash_map<int,vector<int> >& outlist);
		void printMap(hash_map<int,vector<int> >& inlist, hash_map<int,vector<int> >& outlist);

};	

#endif
