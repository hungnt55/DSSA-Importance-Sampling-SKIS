#include <cstdio>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <random>

using namespace std;

int main(int argc, char ** argv)
{
	srand(time(NULL));
	ifstream in(argv[1]);
	unsigned int flag = atoi(argv[3]);
	unsigned long long n,m,i;
	in >> n >> m;
	printf("%lld, %lld\n", n, m);
	vector<unsigned long long> node(n + 1,0);
	vector<unsigned long long> v1,v2,w;
	v1.reserve(m+1);
	v2.reserve(m+1);
	w.reserve(m+1);
	
	unsigned long long t1,t2,pt1,pt2;
	in >> pt1 >> pt2;
	v1.push_back(pt1);
	v2.push_back(pt2);
	w.push_back(1);
	node[pt2]++;
	if (flag == 0) node[pt1]++;
	printf("Reading the graph!\n");
	for (i = 1; i <= m; ++i){
		in >> t1 >> t2;
		node[t2]++;
		if (flag == 0) node[t1]++;

		if(t1 == pt1 && t2 == pt2){
			w[w.size() - 1]++;
			cout << "duplicate " << t1 << " " << t2 << endl;
		} else {
			v1.push_back(t1);
		        v2.push_back(t2);
		        w.push_back(1);
			pt1=t1;
			pt2=t2;
		}
		if (i %100000 == 0)
			printf("%lld\n", i);
	}
	in.close();

	vector<float> rand_n(n+1);
	for (i = 0; i < n+1; ++i){
		rand_n[i] = rand()*1.0/RAND_MAX;
	}

	ofstream out(argv[2]);
	printf("Writing down to file!\n");
	out << n << " ";

	if (flag == 0){
		out << v1.size()*2 << endl;
	}else{
		out << v1.size() << endl;
	}

	for (i = 0; i < v1.size(); ++i){
		out << v1[i] << " " << v2[i] << " " << (float) w[i]*rand_n[v2[i]]/(float)node[v2[i]] << endl;
		if (flag == 0)
			out << v2[i] << " " << v1[i] << " " << (float) w[i]*rand_n[v1[i]]/(float)node[v1[i]] << endl;
	}
	out.close();
}
