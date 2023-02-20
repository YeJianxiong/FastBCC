#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <tuple>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <vector>
#include <tuple>
#include <stdint.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <cilk/cilk_api.h>
#include <cilk/cilk.h>
#include "CycleTimer.h"

//#define DEBUG
//#define PRINTRESULT

typedef unsigned long long DataType;

using namespace std;

template <typename T> class VQueue {
	size_t V;
	size_t curr_;
	size_t next_;
	size_t end_;

public:
	vector<T> data;
	vector<T> value;
	explicit VQueue(size_t V) : data(V), value(V), V(V), curr_(0), next_(0), end_(0){
	}

	inline bool empty() const { return curr_ == next_; }
	inline bool full() const { return end_ == V; }
	inline T &front() { return data[curr_];}
	inline size_t size() const { return end_; }
	inline void pop() { ++curr_; assert(curr_ <= end_);}
	inline void push(const T &val){ data[end_++] = val; assert(end_ <= V);}
	inline void push(const T &dat, const T &val){ data[end_] = dat;value[end_++] = val; assert(end_ <= V);}
	inline void next() {  assert(curr_ == next_);next_ = end_; }//
	inline void clear() { curr_ = next_ = end_ = 0; }
	inline void resize(size_t V_){
		if (V_ > V){ V = V_; data.resize(V); value.resize(V);}
	}
	inline size_t sum(){
		size_t s=0;
		for (DataType i = 0; i < end_; i++){
			s += value[i];
		}
		return s;
	}

	inline typename vector<T>::iterator begin() { return data.begin();}
	inline typename vector<T>::iterator end() { return data.begin() + end_;}
};


vector<vector<DataType>>			Graph;  // the structure of sampled graph
vector<vector<DataType>>			Maps;  // visited states of vertices for all threads
vector<VQueue<DataType>> 			Queues;  // the queue of BFS for all threads
vector<DataType>					Distances;  // for all nodes
vector<DataType>                    reached_Num;  // store the number of source vertex's reachable number for all threads
DataType  							Node_Num = 0, Edge_Num=0, Sample_nodenum=0;
bool 								is_directed = false;
vector<DataType>                    Qset;  // save source vertices
vector<double>                      BCC, SumBCC; // the result of closeness centrality through BFS
vector<pair<DataType, DataType>>    Edge_Origin;  // the edges list of original graph
vector<double>                      Weight;  // the weight(probability) of edges
map<DataType, DataType>             ID2Pos, Pos2ID;  // map the ID of vertex to logical position and the logical position to ID
double                              TIMELIMIT = 20000;  // the time limit is 20000s

DataType num_threads = __cilkrts_get_nworkers();

static inline void ResetMap(DataType threadID)
{
	for (DataType v : Queues[threadID].data) {
		Maps[threadID][v>>5] = 0;
	}
	//Queues[threadID].clear();
}

//set bit at v in Maps
static inline void SetBit(DataType threadID, DataType v)
{
	Maps[threadID][v >> 5] |= (1 << ( v & 31 ));
}

//test bit at v in Maps
static inline int TestBit(DataType threadID, DataType v)
{
	return (Maps[threadID][v >> 5] >> (v & 31) ) & 1;
}

void ReadEdges(string filename, DataType& Qnum)
{
    #ifdef DEBUG
    cerr << "Reading the Edges...\n" ;
    #endif // DEBUG

    //1.fetch edge set from stdin
	FILE *fp = stdin; size_t bufsize=0; char *line = NULL;
	fp = fopen(filename.c_str(),"r");
    int res;
    DataType u, v;

    set<pair<DataType, DataType>> tmpset;
    set<DataType> nodeset;
    while (!feof(fp)){
		res = getline(&line,&bufsize,fp);
		if (res == -1) break;
		if (line[0] == '#') continue;

		res = sscanf(line, "%llu|%llu", &u, &v);
		if ( !res || res == EOF ) {
			continue;
		}

		if(!tmpset.count({u,v})) {
            Edge_Origin.emplace_back(make_pair(u,v));
            tmpset.insert({u,v});
		}
		if(!nodeset.count(u)) {
            nodeset.insert(u);
            //if(Qset.size() < Qnum) Qset.emplace_back(u);
		}
		if(!nodeset.count(v)) {
            nodeset.insert(v);
            //if(Qset.size() < Qnum) Qset.emplace_back(v);
		}
    }
	fclose(fp);
	cin.clear();
	Edge_Num = Edge_Origin.size();
	Node_Num = nodeset.size();

	random_device rd;
    mt19937 g(rd());
    uniform_int_distribution<> dis(0, Node_Num-1);
    set<int> Qtmp;
	while(Qset.size() < Qnum){
        int VID = dis(g);
        if(!Qtmp.count(VID)){
            Qtmp.insert(VID);
            Qset.emplace_back(VID);
        }
	}

	#ifdef DEBUG
    cerr << " Num of nodes " << Node_Num << ", num of edges " << Edge_Num << endl;
    #endif // DEBUG
}

void ReadWeight(string filename, DataType edgenum)
{
    #ifdef DEBUG
    cerr << "Reading the Edge Weight...\n" ;
    #endif // DEBUG

	FILE *fp = stdin; size_t bufsize=0; char *line = NULL;
	fp = fopen(filename.c_str(),"r");
    int res;
    double we;
    DataType num = 0;

    while (!feof(fp) && num <= edgenum){
		res = getline(&line,&bufsize,fp);
		if (res == -1) break;
		if (line[0] == '#') continue;

		res = sscanf(line, "%lf", &we);
		if ( !res || res == EOF ) {
			continue;
		}
        Weight.emplace_back(we);
        num++;
    }
	fclose(fp);
	cin.clear();
}

void SamlpeGraph()
{
    random_device rd;
    mt19937 g(rd());
    uniform_int_distribution<> dis(1, 100);
    Graph.resize(Node_Num);

    //1. sample the edges
    DataType position = 0;
    for(int i=0; i<Edge_Origin.size(); i++){
        if(i>Weight.size()) {cout<<"i>Weight.size()"<<endl; exit(0);}
        double wei = Weight[i];
        double pro = float(dis(g)) / 100;
        if(pro > wei) continue;

        pair<DataType, DataType> edge = Edge_Origin[i];
        auto u = edge.first;
        auto v = edge.second;

        if(!ID2Pos.count(u)){
            ID2Pos[u] = position;
            Pos2ID[position] = u;
            position++;
        }
        if(!ID2Pos.count(v)){
            ID2Pos[v] = position;
            Pos2ID[position] = v;
            position++;
        }
        Graph[ID2Pos[u]].emplace_back(ID2Pos[v]);
		if (!is_directed){
            Graph[ID2Pos[v]].emplace_back(ID2Pos[u]);
        }
    }
    Graph.resize(position);


    //2.sort adjacent lists
	cilk_for (DataType v = 0; v < Graph.size(); v++){
		sort(Graph[v].begin(), Graph[v].end());
	}
	Sample_nodenum = Graph.size();


    //3. Init the graph
    Queues.clear();
    Maps.clear();
    for (DataType t = 0; t < num_threads; t++){
    	Queues.emplace_back(VQueue<DataType>(Sample_nodenum));
    	Maps.emplace_back(vector<DataType>(Sample_nodenum / 32 + 1) );
    	reached_Num.emplace_back(0);
    }
    Distances.resize(Sample_nodenum);
}

//Compute the distance from v to other nodes in networks by BFS
DataType SSSP(DataType v)
{
	DataType threadID = __cilkrts_get_worker_number();

	#ifdef DEBUG
    cout<<"v="<<v<<",";
	cout<<"threadID="<<threadID<<endl;
    #endif // DEBUG

	auto &Q = Queues[threadID];
	int distance = 0;

	Q.clear();
	Q.push(v,0);
	Q.next();
	SetBit(threadID, v);

	while (!Q.empty()) {
		distance++;
		while (!Q.empty()) {
			DataType s = Q.front();
			Q.pop();
			for (DataType w : Graph[s]) {
				if (!TestBit(threadID,w)) {
					Q.push(w, distance);  // Q.next();
					SetBit(threadID,w);
					reached_Num[threadID]++;
				}
			}
		}
		Q.next();
	}

	ResetMap(threadID);

	#ifdef DEBUG
	cout<<"SSSP is complete."<<endl;
    #endif // DEBUG

	return Q.sum();
}

//Compute the BCC
void ComputingBCC(DataType& Q_num)
{
    #ifdef DEBUG
    cout << "ComputingBCC..." << endl;
    #endif // DEBUG

	BCC.resize(Q_num);
	Distances.resize(Q_num);
	vector<DataType> Rnum(Q_num, 0);
	cilk_for (DataType i = 0; i < Q_num; i++){
	    auto VID = Qset[i];
        if(!ID2Pos.count(VID)) {
            Distances[i] = 0;
        }
        else{
            DataType v = ID2Pos[VID];

            Distances[i] = SSSP(v); /// compute distance from source vertex to other vertex

            for (DataType t = 0; t < num_threads; t++){
                Rnum[i] += reached_Num[t];
                reached_Num[t] = 0;
            }
        }

	}

	for (DataType i = 0; i < Q_num; i++){
		if ( Distances[i]==0 || Node_Num<=1) {
            BCC[i] = 0;
		}
		else {
            #ifdef DEBUG
            //cout << "Rnum="<<Rnum[i]<<",Distance="<<Distances[i] << endl;
            #endif // DEBUG
            BCC[i] = Rnum[i] * Rnum[i] * 1.0 / (Node_Num-1) / Distances[i];
		}

		#ifdef DEBUG
		cout<<"BCC["<<Qset[i]<<"]"<<BCC[i]<<endl;
        #endif // DEBUG
	}
	#ifdef DEBUG
    cout << "ComputingBCC complete." << endl;
    #endif // DEBUG
}

void PrintBCC(DataType& Q_num)
{
    cout<<"PrintBCC..."<<endl;
    for(DataType i=0; i<Q_num; i++){
        DataType VID = Qset[i];
        cout<<"BCC[ "<<VID<<" ]= "<<SumBCC[i]<<endl;
    }
    cout<<endl;
}

void PrintGraph()
{
    cout<<"PrintGraph..."<<endl;
    for(int i = 0; i < Graph.size(); i++){
        cout<<Pos2ID[i]<<": ";
        for(auto VPos: Graph[i]){
            DataType v = Pos2ID[VPos];
            cout<<v<<" ";
        }
        cout<<endl;
    }
}

void Clear()
{
    Graph.clear();
    ID2Pos.clear();
    Pos2ID.clear();
    Distances.clear();
    BCC.clear();
}

int main(int argc, char *argv[])
{
	if (argc < 5) {
        cout << "need more args..." << endl;
		cout << "executation filename directed Q_filepath" << endl;
        return 0;
    }

    string dataset_file = argv[1];
	is_directed = atoi(argv[2]);

	//Build(argv[1]);
    //PrintGraph();

    DataType Q_num = atoi(argv[3]);
    ReadEdges(dataset_file, Q_num);
    Q_num = Qset.size();

    string weight_file = argv[4];
    ReadWeight(weight_file, Edge_Num);

    if(argc >= 6) num_threads = atoi(argv[5]);
	__cilkrts_set_param("nworkers",to_string(num_threads).data());

	int sample_num = 1024;
	if(argc == 7) sample_num = atoi(argv[6]);


	#ifdef DEBUG
    cout<<"Dataset is "<<argv[1]<<endl;
    cout<<"The graph is directed or not (1-direted/0-undirected):"<<is_directed<<endl;
    cout<<"|Q| = "<<Q_num<<endl;
    cout << "weigth_file = " << weight_file << endl;
    cout << "sample_num = " << sample_num << endl;
    cout << "NumThread = " << num_threads << endl;
    #endif // DEBUG

    double sample_times, BCC_times, start, Total_times;
    double time_limit = CycleTimer::currentSeconds() + TIMELIMIT;
    bool isTLE = false;

    sample_times = 0;
    BCC_times = 0;
    SumBCC.resize(Q_num);
	for(int ss=0; ss<sample_num; ss++){
        if(CycleTimer::currentSeconds() > time_limit){
            isTLE = true;
            break;
        }
        #ifdef DEBUG
        cout<<"The "<<ss<<"th sample..."<<endl;
        #endif // DEBUG

        start = CycleTimer::currentSeconds();
        SamlpeGraph();
        sample_times += CycleTimer::currentSeconds() - start;

        //PrintGraph();

        start = CycleTimer::currentSeconds();
        ComputingBCC(Q_num);
        if(ss==0) {
            for(int i=0; i<Q_num; i++){
                SumBCC[i] = BCC[i];
            }
        }
        else{
            for(int i=0; i<Q_num; i++){
                SumBCC[i] += BCC[i];
            }
        }
        BCC_times += CycleTimer::currentSeconds() - start;
        Clear();

        #ifdef DEBUG
        cout<<"The "<<ss<<"th sample is done.\n\n"<<endl;
        #endif // DEBUG
	}

	start = CycleTimer::currentSeconds();
    for(int i=0; i<Q_num; i++){
        SumBCC[i] /= sample_num;
    }
	BCC_times += CycleTimer::currentSeconds() - start;
    Total_times = sample_times + BCC_times;

    #ifdef DEBUG
    cout << "Sample Time: " << sample_times << "s" << endl;
    cout << "baselineBCC Time: " << BCC_times << "s" << endl;
    cout << "Total time: " << Total_times << "s" << endl;
    #endif // DEBUG

    #ifdef PRINTRESULT
    PrintBCC(Q_num);
    #endif // DEBUG

    FILE* fp;
    fp = freopen("result-FastBCC.txt","a", stdout);
    if(isTLE){
        cout<<"["<<dataset_file<<"] "<<"TLE"<<" "<<"TLE"<<" "<<"TLE"<<" "<<Q_num<<" "<<sample_num<<" "<<num_threads<<" FastBCC"<<endl;
    }
    else{
        cout<<"["<<dataset_file<<"] "<<Total_times<<" "<<BCC_times<<" "<<sample_times<<" "<<Q_num<<" "<<sample_num<<" "<<num_threads<<" FastBCC"<<endl;
    }
    fclose(fp);

	return 0;
}
/**
Readme

./FastBCC ../test_queries/data/test.txt 1 5 ./weight/sample_test.txt
./FastBCC ../test_queries/data/test.txt 1 5 ./weight/sample_test.txt 1
./FastBCC ../test_queries/data/test.txt 1 10 ./weight/sample_test.txt 8 1024

input:
[executaion file] [datasets file] [is_direccted] [number of source vertices] [file of weight] ([number of thread]) ([sample number])

output:
[datasets file] [total time] [computation time] [sample time] [number of source vertices] [sample number] [thread number] [algorithm name]

*/

