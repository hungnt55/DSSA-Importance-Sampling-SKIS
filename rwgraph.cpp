#include "rwgraph.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <cstdlib>
#include <unistd.h>
#include <sstream>

using namespace std;

double Graph::getTotalProb()
{
	return totalProb;
}

const vector<int> & Graph::operator [] (int u) const
{
	return adjList[u];
}


const vector<int> & Graph::operator [] (int u)
{
	return adjList[u];
}


const vector<UI> & Graph::getWeight (int u) const
{
        return weights[u];
}

const vector<UI> & Graph::getWeight (int u)
{
        return weights[u];
}

/*
* get degree of node u
*/
int Graph::getDegree(int u) const
{
	return adjList[u].size();
}

/*
* get the number of nodes
*/
int Graph::getSize() const
{
	return numNodes;
}

/*
* get the number of edges
*/
int Graph::getEdge() const
{
	return numEdges;
}

void Graph::readGraphLT(const char * filename)
{
	FILE * pFile;
        pFile = fopen(filename, "rb");
        fread(&numNodes, sizeof(int), 1, pFile);
        fread(&numEdges, sizeof(long long), 1, pFile);
        node_deg=vector<int>(numNodes + 1);
        fread(&node_deg[1], sizeof(int), numNodes, pFile);
        
        vector<int> a;
        vector<UI> b;
        adjList.push_back(a);
        weights.push_back(b);

        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<int> tmp(node_deg[i]);
                fread(&tmp[0], sizeof(int), node_deg[i], pFile);
                adjList.push_back(tmp);
        }
                
        totalProb = 0;
        nodeProb = vector<UI>(numNodes+1,0);
        nonProb = vector<float>(numNodes+1,0);
        vector<double> nodeProbT(numNodes+1,0);
        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<float> tmp(node_deg[i] + 1, 0);
                vector<UI> tmp1(node_deg[i] + 1, 0);
                fread(&tmp[1], sizeof(float), node_deg[i], pFile);

//		cout << "Neighbors: ";
                for(int j = 1; j < node_deg[i]+1; ++j){
                        nonProb[i] += tmp[j];
			tmp1[j] = nonProb[i]*UI_MAX;
			if (nonProb[i] >= 1){
				tmp1[j] = UI_MAX;
			}
//			cout << tmp1[j] << " ";
                }
//		cout << endl;
                totalProb += nonProb[i];
		nonProb[i] = 1 - nonProb[i];
                nodeProbT[i] = totalProb;

                weights.push_back(tmp1);
        }
        for (unsigned int i = 1; i <= numNodes; ++i){
                nodeProb[i] = (nodeProbT[i]/totalProb)*UI_MAX;
        }

        vector<double>().swap(nodeProbT);
}

/*
* read input graph for IC model
*/
void Graph::readGraphIC(const char* filename)
{
    	FILE * pFile;
    	pFile = fopen(filename, "rb");
    	fread(&numNodes, sizeof(int), 1, pFile);
    	fread(&numEdges, sizeof(long long), 1, pFile);
    	node_deg=vector<int>(numNodes + 1);
    	fread(&node_deg[1], sizeof(int), numNodes, pFile);
        
	vector<int> a;
	vector<UI> b;
    	adjList.push_back(a);
    	weights.push_back(b);
	
        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<int> tmp(node_deg[i]);
                fread(&tmp[0], sizeof(int), node_deg[i], pFile);
                adjList.push_back(tmp);
        }

	totalProb = 0;
	nodeProb = vector<UI>(numNodes+1,0);
	nonProb = vector<float>(numNodes+1,1);
	vector<double> nodeProbT(numNodes+1,0);
        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<float> tmp(node_deg[i] + 1, 0);
		vector<UI> tmp1(node_deg[i] + 1, 0);
                fread(&tmp[1], sizeof(float), node_deg[i], pFile);

                for(int j = 1;j < node_deg[i] + 1; ++j){
                        tmp1[j] = (tmp[j]*nonProb[i])*UI_MAX + tmp1[j-1];
			nonProb[i] *= (1-tmp[j]);
                }
		totalProb += 1-nonProb[i];
		nodeProbT[i] = totalProb;

		if (tmp1[node_deg[i]] <= 0)
                        tmp1[node_deg[i]] = UI_MAX;
		
                weights.push_back(tmp1);
        }
	for (unsigned int i = 1; i <= numNodes; ++i){
		nodeProb[i] = (nodeProbT[i]/totalProb)*UI_MAX;
	}
	
	vector<double>().swap(nodeProbT);
}

void Graph::writeToFile(const char * filename)
{
}

// choose a random edge in LT model based on linear search
inline int HyperGraph::randIndex_lin(const vector<UI> &w, unsigned int si, sfmt_t &sfmtSeed)
{
        UI ranNum = sfmt_genrand_uint32(&sfmtSeed);
        if (si <= 1 || ranNum > w[si - 1])
                return -1;

        for (unsigned int i = 1; i < si; ++i){
                if (ranNum <= w[i])
                        return i;
        }
        return -1;
}

// binary search is used
inline int HyperGraph::randIndex_bin(const vector<UI> &w, unsigned int si, sfmt_t &sfmtSeed, int ben, int end, bool normalized)
{
        UI ran = sfmt_genrand_uint32(&sfmtSeed);
	// cout << "Rand: " << ran << endl << endl << endl;
        if (normalized){
                ran *= w[end-1]/(double)UI_MAX;
        } else {
                ran = ran*((UI_MAX-w[ben])/(double)UI_MAX) + w[ben];
        }
        if (si <= 1 || ran > w[end - 1] || ben == end-1)
                return -1;
        int left = ben;
        int right = end - 1;
        int prob = 0;
        for (unsigned int i = 0; i < si; ++i){
                prob = (left + right)/2;
		if (w[prob] <= ran){
                        left = prob + 1;
                        continue;
                }
                if (w[prob - 1] > ran){
                        right = prob - 1;
                        continue;
                }
                break;
        }
        return prob;
}

HyperGraph::HyperGraph(Graph & g, unsigned int n, unsigned int m)
{
	sfmt_init_gen_rand(&sfmtSeed, rand());
	node_edge = vector<vector<int> >(n+1);
	maxDegree = 0;
	numNodes = n;
	numEdges = m;
	curEdge=0;
	nonProb = g.nonProb;
	gamma = g.totalProb;
}

vector<float> & HyperGraph::getNonProb()
{
	return nonProb;
}

double HyperGraph::getGamma()
{
	return gamma;
}

void HyperGraph::updateDeg(){
	unsigned int num=edge_node.size();
	for (unsigned int i = curEdge; i < num; ++i){
		unsigned int num2 = edge_node[i].size();
		for (unsigned int j=0;j<num2;++j){
			node_edge[edge_node[i][j]].push_back(i);
		}
	}
	curEdge = edge_node.size();
}

void HyperGraph::updateEdge(){
	curEdge = edge_node.size();
}

/*
* Add a hyperedge into the hypergraph
*/
void HyperGraph::addEdge(vector<int> & edge)
{
	edge_node.push_back(edge);
	unsigned int ind = edge_node.size() - 1;
	for (unsigned int i = 0; i < edge.size(); ++i)
		node_edge[edge[i]].push_back(ind);
}

/*
* Add a hyperedge into the hypergraph while keeping track of the node with max degree
*/
void HyperGraph::addEdgeD(vector<int> & edge)
{
        edge_node.push_back(edge);
        int ind = edge_node.size() - 1;
        for (unsigned int i = 0; i < edge.size(); ++i){
                node_edge[edge[i]].push_back(ind);
		if (node_edge[edge[i]].size() > maxDegree)
			maxDegree = node_edge[edge[i]].size();
	}
}

/*
* get an edge from the hypergraph
*/
const vector<int> & HyperGraph::getEdge(int e) const{
	return edge_node[e];
}

const vector<int> & HyperGraph::getEdge(int e){
	return edge_node[e];
}

/*
* get the list of hyperedges incident to node n
*/
const vector<int> & HyperGraph::getNode(int n) const{
	return node_edge[n];
}

const vector<int> & HyperGraph::getNode(int n){
	return node_edge[n];
}

/*
* get the number of hyperedges
*/
int HyperGraph::getNumEdge() const
{
        return edge_node.size();
}

/*
* get the maximum degree
*/
int HyperGraph::getMaxDegree()
{
	return maxDegree;
}

/*
* remove all the hyperedges
*/
void HyperGraph::clearEdges()
{
	edge_node.clear();
	node_edge.clear();
	cout << "clear edges!" << endl;
	maxDegree = 0;
}

void HyperGraph::pushHyperedges(vector<vector<int> > & hyperedges){
	for (unsigned int i = 0; i < hyperedges.size(); ++i){
		edge_node.push_back(hyperedges[i]);
	}
}

// polling process under IC model which also checks whether it intersects with a seed set 
bool HyperGraph::pollingIC(Graph &g, vector<bool> & link, vector<bool> &visit, vector<int> &visit_mark, sfmt_t & sfmtSeed1)
{
        int selectedPos;
	unsigned int cur = randIndex_bin(g.nodeProb,g.numNodes+1,sfmtSeed1,0,g.numNodes+1,true);
	if (link[cur]){
		return true;
	}
	bool returnedValue = false;
        unsigned int num_marked=0;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;

	const vector<UI> &w=g.getWeight(cur);
        int wsize = w.size();
        const vector<int> &neigh = g[cur];
        selectedPos = randIndex_bin(w,wsize,sfmtSeed1,0,wsize,true);
        while (selectedPos > -1){
                if (!visit[neigh[selectedPos-1]]){
			if (link[neigh[selectedPos-1]]){
				returnedValue = true;
				break;
			}

                        visit[neigh[selectedPos-1]] = true;
                        visit_mark[num_marked]=neigh[selectedPos-1];
                        num_marked++;
                 }
                 selectedPos = randIndex_bin(w,wsize,sfmtSeed1,selectedPos,wsize,false);
        }

        unsigned int curPos=1;
        while(curPos < num_marked && !returnedValue){
                cur = visit_mark[curPos];
                curPos++;
                const vector<UI> &w=g.getWeight(cur);
                int wsize = w.size();
                const vector<int> &neigh = g[cur];
                selectedPos = randIndex_bin(w,wsize,sfmtSeed1,0,wsize,false);
                while (selectedPos > -1){
			if (link[neigh[selectedPos-1]]){
				returnedValue = true;
				break;
			}
                        if (!visit[neigh[selectedPos-1]]){
                                visit[neigh[selectedPos-1]] = true;
                                visit_mark[num_marked]=neigh[selectedPos-1];
                                num_marked++;
                        }
                        selectedPos = randIndex_bin(w,wsize,sfmtSeed1,selectedPos,wsize,false);
                }
        }
//	cout << "Sample size: " << num_marked << endl;

        for(unsigned int i = 0; i < num_marked;++i){
                visit[visit_mark[i]]=false;
        }
	
	return returnedValue;
}


// Pure generating samples under the IC model
void HyperGraph::pollingIC1(Graph &g, vector<bool> &visit, vector<int> &visit_mark, sfmt_t & sfmtSeed1)
{
        int selectedPos;
        unsigned int cur = randIndex_bin(g.nodeProb,g.numNodes+1,sfmtSeed1,0,g.numNodes+1,true);
	// cout << "Source: " << cur << endl;
        unsigned int num_marked = 0;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;

        const vector<UI> &w=g.getWeight(cur);
        int wsize = w.size();
        const vector<int> &neigh = g[cur];
        selectedPos = randIndex_bin(w,wsize,sfmtSeed1,0,wsize,true);
        while (selectedPos > -1){
                if (!visit[neigh[selectedPos-1]]){
                        visit[neigh[selectedPos-1]] = true;
                        visit_mark[num_marked]=neigh[selectedPos-1];
                        num_marked++;
                 }
                 selectedPos = randIndex_bin(w,wsize,sfmtSeed1,selectedPos,wsize,false);
        }

        unsigned int curPos=1;
        while(curPos < num_marked){
                cur = visit_mark[curPos];
                curPos++;
                const vector<UI> &w=g.getWeight(cur);
                int wsize = w.size();
                const vector<int> &neigh = g[cur];
                selectedPos = randIndex_bin(w,wsize,sfmtSeed1,0,wsize,false);
                while (selectedPos > -1){
                        if (!visit[neigh[selectedPos-1]]){
                                visit[neigh[selectedPos-1]] = true;
                                visit_mark[num_marked]=neigh[selectedPos-1];
                                num_marked++;
                        }
                        selectedPos = randIndex_bin(w,wsize,sfmtSeed1,selectedPos,wsize,false);
                }
        }
//	cout << "Sample size: " << num_marked << endl;
//	cout << "Sizes: " << num_marked << endl;
//	cout << "RR set: ";
        for(unsigned int i = 0; i < num_marked;++i){
//		cout << visit_mark[i] << " ";
                visit[visit_mark[i]]=false;
        }
//	cout << endl;
        edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
}

bool HyperGraph::pollingIC2(Graph &g, vector<bool> & link, vector<bool> &visit, vector<int> &visit_mark, sfmt_t & sfmtSeed1)
{
        int selectedPos;
        unsigned int cur = randIndex_bin(g.nodeProb,g.numNodes+1,sfmtSeed1,0,g.numNodes+1,true);
//	cout << "Source: " << cur << endl;
        bool returnedValue = false;

        if (link[cur]){
                returnedValue = true;
        }

        unsigned int num_marked=0;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;
        const vector<UI> &w=g.getWeight(cur);
        int wsize = w.size();
        const vector<int> &neigh = g[cur];
        selectedPos = randIndex_bin(w,wsize,sfmtSeed1,0,wsize,true);

        while (selectedPos > -1){
                if (!visit[neigh[selectedPos-1]]){
                        if (link[neigh[selectedPos-1]]){
                                returnedValue = true;
                        }
                        visit[neigh[selectedPos-1]] = true;
                        visit_mark[num_marked]=neigh[selectedPos-1];
                        num_marked++;
                 }
                 selectedPos = randIndex_bin(w,wsize,sfmtSeed1,selectedPos,wsize,false);
        }

	unsigned int curPos=1;
        while(curPos < num_marked){
                cur = visit_mark[curPos];
                curPos++;
                const vector<UI> &w=g.getWeight(cur);
                int wsize = w.size();
                const vector<int> &neigh = g[cur];
                selectedPos = randIndex_bin(w,wsize,sfmtSeed1,0,wsize,false);

                while (selectedPos > -1){
                        if (link[neigh[selectedPos-1]]){
                                returnedValue = true;
                        }
                        if (!visit[neigh[selectedPos-1]]){
                                visit[neigh[selectedPos-1]] = true;
                                visit_mark[num_marked]=neigh[selectedPos-1];
                                num_marked++;
                        }
                        selectedPos = randIndex_bin(w,wsize,sfmtSeed1,selectedPos,wsize,false);
                }
        }
//	cout << "Sample size: " << num_marked << endl;
//	cout << "Sizes: " << num_marked << endl;
//	cout << "RR set: ";
        for(unsigned int i = 0; i < num_marked;++i){
                visit[visit_mark[i]]=false;
//		cout << visit_mark[i] << " ";
        }
//	cout << endl;
	edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));

        return returnedValue;
}

bool HyperGraph::pollingLT(Graph &g, std::vector<bool> & link, std::vector<bool> &visit, std::vector<int> &visit_mark, sfmt_t & sfmtSeed1)
{
	unsigned int i;
        bool t = false;
        unsigned int gSize = g.getSize();
        unsigned int cur = randIndex_bin(g.nodeProb,g.numNodes+1,sfmtSeed1,0,g.numNodes+1,true);
        unsigned int num_marked = 0;
	
	if (link[cur]){
                return true;
        }

        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;

	int ind = randIndex_bin(g.weights[cur],g.node_deg[cur]+1,sfmtSeed1,0,g.node_deg[cur]+1,true);

        cur = g.adjList[cur][ind-1];

        for (i = 0; i < gSize; ++i){
                if (link[cur]){
                        t=true;
                        break;
                }
                if (visit[cur] == true) break;
                visit[cur] = true;
                visit_mark[num_marked] = cur;
                num_marked++;

                ind = randIndex_bin(g.weights[cur],g.node_deg[cur]+1,sfmtSeed1,0,g.node_deg[cur]+1,false);

                if (ind == -1)
                        break;

                cur = g.adjList[cur][ind-1];
        }
        //cout << "Sample size: " << num_marked << endl;
        for (i = 0; i < num_marked; ++i){
                visit[visit_mark[i]]=false;
        }
        return t;
}

void HyperGraph::pollingLT1(Graph &g, std::vector<bool> &visit, std::vector<int> &visit_mark, sfmt_t & sfmtSeed1)
{
	unsigned int i;
        unsigned int gSize = g.getSize();
        unsigned int cur = randIndex_bin(g.nodeProb,g.numNodes+1,sfmtSeed1,0,g.numNodes+1,true);
	cout << cur << ":" << g.node_deg[cur] << " ";
        unsigned int num_marked = 0;

        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;

        int ind = randIndex_bin(g.weights[cur],g.node_deg[cur]+1,sfmtSeed1,0,g.node_deg[cur]+1,true);

        cur = g[cur][ind-1];
//	cout << ind << ":" << cur << " ";

        for (i = 0; i < gSize; ++i){
                if (visit[cur] == true) break;
                visit[cur] = true;
                visit_mark[num_marked] = cur;
                num_marked++;
                const vector<int> &neigh = g[cur];

                ind = randIndex_bin(g.weights[cur],g.node_deg[cur]+1,sfmtSeed1,0,g.node_deg[cur]+1,false);

                if (ind == -1)
                        break;

                cur = neigh[ind-1];
        }
        //cout << "Sample size: " << num_marked << endl;
        edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
        for (i = 0; i < num_marked; ++i){
//		cout << visit_mark[i] << " ";
                visit[visit_mark[i]]=false;
        }
//	cout << endl;
}

bool HyperGraph::pollingLT2(Graph &g, std::vector<bool> & link, std::vector<bool> &visit, std::vector<int> &visit_mark, sfmt_t & sfmtSeed1)
{
	unsigned int i;
        bool t = false;
        unsigned int gSize = g.getSize();
        unsigned int cur = randIndex_bin(g.nodeProb,g.numNodes+1,sfmtSeed1,0,g.numNodes+1,true);
        unsigned int num_marked = 0;

	if (link[cur]){
                t = true;
        }

        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;

        int ind = randIndex_bin(g.weights[cur],g.node_deg[cur]+1,sfmtSeed1,0,g.node_deg[cur]+1,true);

        cur = g[cur][ind-1];

        for (i = 0; i < gSize; ++i){
                if (visit[cur] == true) break;
                visit[cur] = true;
                visit_mark[num_marked] = cur;
                num_marked++;
                if (link[cur])
                        t=true;
                
		ind = randIndex_bin(g.weights[cur],g.node_deg[cur]+1,sfmtSeed1,0,g.node_deg[cur]+1,false);

                if (ind == -1)
                        break;

                cur = g.adjList[cur][ind-1];
        }
        //cout << "Sample size: " << num_marked << endl;

        edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
        for (i = 0; i < num_marked; ++i){
                visit[visit_mark[i]]=false;
        }
        return t;
}

/*
* convert from an integer to a string
*/
string intToStr(int i) {
        stringstream ss;
        ss << i;
        return ss.str();
}

/*
* convert from a strong to an integer
*/
unsigned int strToInt(string s) {
        unsigned int i;
        istringstream myStream(s);

        if (myStream>>i) {
                return i;
        } else {
                cout << "String " << s << " is not a number." << endl;
                return atoi(s.c_str());
        }
        return i;
}
