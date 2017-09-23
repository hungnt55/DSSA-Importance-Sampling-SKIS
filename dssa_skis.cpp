#include "option.h"
#include "hypergraph.hpp"
#include "getMem.hpp"
#include "sfmt/SFMT.h"
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

vector<sfmt_t> sfmtSeed;
bool model = false;

bool calculateInfluence(Graph & g, HyperGraph & hg, vector<int> & seeds, int t, double expected, double epsilon, float delta, long long int maxSamples, double nonProb, double gamma, int iter){
	//cout << "Calculate Influence" << endl;
	long long  counter = 0;
	int n = g.getSize();
	unsigned int k = seeds.size();
	double f = log2(n*(log(2/delta)+lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))/(k*log(6*log2(n)/delta)));
	//cout << f << endl;
	double lambda1 = 1+(2+1/3)*(1+1/2)*(log(3/delta)+log(f))*4;
	vector<bool> link(n+1,false);
	double degree=0;
	for (unsigned int i = 0; i < k;++i){
		link[seeds[i]] = true;
	}
	vector<bool> maxSeed(t, false);

	omp_set_num_threads(t);
	#pragma omp parallel
	{
		vector<bool> visit(n+1,false);
		vector<int> visit_mark(n,0);
		int id = omp_get_thread_num();
		while (counter < maxSamples){
			if (model){
				maxSeed[id]=hg.pollingIC2(g,link,visit,visit_mark,sfmtSeed[id]);
			} else {
				maxSeed[id]=hg.pollingLT2(g,link,visit,visit_mark,sfmtSeed[id]);
			}
                        #pragma omp critical
                        {
                                counter ++;
                                if (maxSeed[id]){
                                        degree ++;
                                }
                        }
                }
	}
	degree += nonProb*counter/gamma;
	//cout << degree << " " << counter << " " << expected << " " << (degree*gamma/counter) << endl;
	
	if (degree >= lambda1){
                double epsilon_1 = (expected)/(degree*gamma/counter) - 1;
                //cout << "Epsilon_1 = " << epsilon_1 << endl;
                double epsilon_2 = sqrt(n*(1+1/2)/(degree*gamma*pow(2,iter-1)/counter));
                //cout << "Epsilon_2 = " << epsilon_2 << endl;
                double epsilon_3 = sqrt(n*(1+1/2)*(1-1/exp(1)-epsilon)/((1+epsilon/3)*degree*gamma*pow(2,iter-1)/counter));
                //cout << "Epsilon_3 = " << epsilon_3 << endl;
                //cout << "Epsilon_t = " << (epsilon_1 + epsilon_2 + epsilon_1*epsilon_2)*(1-1/exp(1)-epsilon) + epsilon_3*(1 - 1/exp(1)) << endl;

                if ((epsilon_1 + epsilon_2 + epsilon_1*epsilon_2)*(1-1/exp(1)-epsilon) + epsilon_3*(1 - 1/exp(1)) <= epsilon){
                        return true;
                }
        }

        hg.updateDeg();
        return false;
}

int main(int argc, char ** argv)
{
	srand(time(NULL));
	
	OptionParser op(argc, argv);
	if (!op.validCheck()){
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}
	char * inFile = op.getPara("-i");
	if (inFile == NULL){
		inFile = (char*)"network.bin";
	}

	char * outFile = op.getPara("-o");
	if (outFile == NULL){
		outFile = (char*)"network.seeds";
	}

        char * mod = op.getPara("-m");
        if (mod != NULL && strcmp(mod,"IC") == 0)
                model = true;

	Graph g;
	if (model){
		g.readGraphIC(inFile);
	} else {
		g.readGraphLT(inFile);
	}

	int n = g.getSize();
	int m = g.getEdge();

	char * tmp = op.getPara("-epsilon");
	float epsilon = 0.1;
	if (tmp != NULL){
		epsilon = atof(tmp);
	}

	float delta = 1.0/n;
	tmp = op.getPara("-delta");
        if (tmp != NULL){
                delta = atof(tmp);
        }

	float k = 100;
	
	tmp = op.getPara("-k");
	if (tmp != NULL){
		k = atof(tmp);
	}	

	int t = 1;
	tmp = op.getPara("-t");
	if (tmp != NULL){
		t = atoi(tmp);
	}

	sfmtSeed = vector<sfmt_t>(t+1);
        for (int i = 0; i <= t; ++i){
                sfmt_init_gen_rand(&sfmtSeed[i], rand());
        }
	
	HyperGraph hg(g,n,m);
	vector<double> degree(k+1,0);

	vector<float> nonProba = hg.getNonProb();

        double f = (log(6/delta)+lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))*n/(k*log(6*log2(n)/delta));

        double lambda = (2+2/3)*log(3*log2(f)/delta);

        long long int totalSamples = (long long int)lambda;
        // cout << f << " " << lambda << " " << totalSamples << endl;

	addHyperedge(g,hg,t,totalSamples,sfmtSeed,model);
	double nmax = (2+2*epsilon/3)*(lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1) + log(6/delta))*n/(epsilon*epsilon*k);

	vector<int> seeds;
	int iter = 1;

	double ob;
	float nonProb = 0;
	
	clock_t start = clock();

	while (totalSamples < nmax){
		seeds.clear();
		totalSamples = hg.getNumEdge();
		//cout << "Find seed set: " << totalSamples << endl << endl << endl;
		buildSeedSet(hg,seeds,n,k,degree,ob);
		nonProb = 0;
		for (int j = 0; j < k; ++j){
			nonProb += nonProba[seeds[j]];
		}
		
		if (calculateInfluence(g,hg,seeds,t,(double)degree[k]*hg.getGamma()/(double)hg.getNumEdge(),epsilon,delta,totalSamples, nonProb, hg.getGamma(),iter)){
                      	break;
                }
		iter++;
	}
	cout << "Total samplings: " << totalSamples << endl;
	cout << "Seed Nodes: ";
	ofstream out(outFile);
	for (unsigned int i = 0; i < seeds.size(); ++i){
		cout << seeds[i] << " ";
		out << seeds[i] << endl;
	}
	out.close();
	cout << endl << endl;
	printf("Influence: %0.2lf\n",(double)degree[k]*hg.getGamma()/(double)totalSamples);
	cout << "Time: " << (float)(clock()-start)/CLOCKS_PER_SEC << "s" << endl;
	cout << "Memory: " << getMemValue()/1024.0 << endl;
}
