#ifndef TIC_H
#define TIC_H

#include "global.h"
#include <vector>
#include <map>
#include <set>

using namespace std;

typedef std::pair<double, int> di;

struct Edge 
{
	int u,v;
	vector<double> w;
};

struct UI
{
    int uid;
    int kid;
    bool operator<(const UI&ui)const{
        if(uid<ui.uid)return true;
        if(uid==ui.uid &&kid<ui.kid)return true;
        return false;
    }
};


class TIC
{
private:
	static int n_;
	static int m_;
	static int Z_;
    static int K_;

    static map<int, int> cvt_;  // ID converter
    static map<int, int> cvt2_;  // ID converter
    static map<int, double> cost_;
    static map<UI,int> ui_;
    static map<int,vector<int> > ui2_;
    static map<int, vector<double> > kt_; 
//	static vector<int> indegree_;
//	static vector<int> outdegree_;
	static vector<int> inindex_;
	static vector<int> outindex_;
	static vector<Edge> inedges_;
	static vector<Edge> outedges_;
	static vector<double> outcache_;
	static vector<double> incache_;
	static vector<bool> outcflag_;  // is cached ?
	static vector<bool> incflag_;

    static void init();

	static void qsort_inedges(int h, int t);
	static void qsort_outedges(int h, int t);
    static void RemoveMultiPaths();
    static void IndexNodesOnEdges();

public:
    static double GetCost(int v);
    static int GetInterest(int,int);
    static vector<int> GetInterest(int);
    static vector<double> GetTopic(int);
    static int GetName(int v);
    static int GetId(int v);
	static int	GetN();
	static int	GetM();
	static int	GetZ();
    static int  GetK();

//	static int	GetOutDegree(int vertex);
//	static int	GetInDegree(int vertex);
	static int	GetOutNeighbor(int vertex);
	static int	GetInNeighbor(int vertex);
    static double* GetOutNeighborBound(int vertex);
    static double* GetInNeighborBound(int vertex);
	static Edge	GetOutEdge(int vertex, int idx);
	static Edge	GetInEdge(int vertex, int idx);
	static double	GetOutEdgeWeight(int vertex, int idx, const vector<double>&);
	static double	GetInEdgeWeight(int vertex, int idx, const vector<double>&);

	static void Build2TICFromFile(char *);
    static void BuildCost(char *);
    static void BuildKeywordTopic(char *);
    static void BuildUserInterest(char *);
    static void Print();
    static void SaveConverter(char *);
    static void Stats();

    static double MaxPProb(vector<double>&);
    static double PProb(const vector<double>&, const vector<double>&);
    static vector<di> Dijkstra(int u, const double theta);
    static vector<di> ReverseDijkstra(int u, const double theta);
    static vector<di> DijkstraWithStopIds(int u, const double theta, const set<int>&);
    static vector<di> DijkstraOffline(int u, const double theta);


    static double Spread(vector<int>, vector<double>);
    static double WeightedScore(vector<int>, vector<double>,vector<int>);
    static double Convert2PMIA(vector<double>, char *);  /* convert into topic-aware PMIA format, return the convert time cost.*/
    static void Convert2NoTopicPMIA(char *);  /* convert into topic-aware PMIA format, return the convert time cost.*/
    static void resetCache();
};

#endif
