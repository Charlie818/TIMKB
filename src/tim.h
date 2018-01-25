#ifndef TIM_H
#define TIM_H

#include <set>
#include <vector>

#include "tic.h"

#define BOUND      0
#define NOBOUND    1
#define INEDGE     2
#define RAWBOUND   3

using namespace std;
typedef std::pair<int,double> tuple;
typedef struct{
	vector<int> keys;
	int n;
}Key;
typedef struct Inverted_List
{
    int kid;
    double idf;
    vector<tuple> list;
}IL;
// nb[kid]->(u,max,ub)
typedef struct NeighborHood_Bound
{
    int uid;
    double max_ub;
    vector<double> ub;
}NB;
typedef struct
{
    vector<double> item;
    double B;
    vector<int> keywords;
} Q;							//query

typedef struct
{
    vector<int> ids;
    vector<double> minflu;
    vector<double> cost;
} S;							//S set

typedef struct
{
    double B;
    double min_ratio;
    double max_ratio;
    double spread;
} BNode;    


#define INITIAL 0
#define ESTIMATED 1
#define COMPUTED 2
typedef struct
{
    int id;
    int name;
    double benefit;
    int status;  // INITIAL,ESTIMATED, COMPUTED
    double cost;
    int round;
    double influ;
    vector<double> vec;
} HNode;

typedef struct
{
    int id;
    double infl;
    vector<double> hist;  // theta / theta2_
    int size;
} LNode;

// Lower Bound
typedef struct
{
//    vector<vector<int> > lbouts;
    vector<int> seeds;
    vector<double> item;
    vector<double>cost;
    vector<double> influ;
} LBSample;

//Upper Bound
typedef struct
{
    vector<double> item;
    vector<BNode> nodes;
} UBSample;

class TIM 
{
private:
    static double GetIDF(int);
    static double GetIDF(vector<int>);
    static void CursorInit(Q );
    
    static double GetKUWeights(int,vector<int> );
  
    
    static vector<vector<NB> > nb_;
    static vector<NB> nb2_;
    static HNode pop_H();
    static void push_H(HNode);
    static di BestFirstWithoutCost(Q);
    static di BestFirstWithCost(Q);
    static di NeighborhoodWithCost(Q q);
    static di GreedyWithCost(Q);
    static void InsertCandidates(vector<di>);
    static void InsertCandidates2(Q q);
    static void InsertCandidatesWithCost(Q);
    static void InsertCandidatesWithoutCost(vector<di>,Q);
    static void InsertCandidatesWithCostGreedy(vector<di> );
    static double EstMarginUB(HNode,Q);
    static double EstMarginUB2(HNode seed,Q q);
    static double CalcMargin(HNode,Q); 
    static vector<HNode> H;
    // needed for all methods
    static int n_;
    static double *ap_;
    static Q q_;
    static int round_;
    static bool *used_;

    // not needed for every method
    // Greedy
    static vector<HNode> H_;
    static vector<vector<di> > iees_cache_;
    static double* tap_;  // temporary ap
//    static vector<int> pre_;
//    static vector<set<int> > In_;
    // List
    static vector<di> L_;  // \List  <influ, number_of_iee>
    static int *size_;  /* size of offline uset */
//    static double *infl2_;  /* sqrt(theta_) */
    static vector<double *>hist_;  /* histogram */
    static double* bars_;  /* theta bars */
    static int cursorL_;  // cursor on \List
    // Dijkstra
    static double *dist_;
    static int *seen_;
    static int *seen_idx_;
    static int *children_;
    static int *pred_;  /* predecessor id */
    static bool *is_border_;  /* if on the border, upper bounded needed */
    /*in-edge*/
    static vector<vector<double> > max_topic_; 
//    static vector<vector<di> > USet_;
//    static vector<vector<double> > max_topic_;
//    static double *max_in_edge_;  // \List  <influ, number_of_iee>

//    static vector<vector<di> > Iee_;  // exact calculated out arborescence
    // Lower
    static double *lap_;
    static vector<LBSample> lbsamples_;
    // upper
    static vector<UBSample> ubsamples_;

    static int Dijkstra(int, double);
    static void UpdateHeap(int);
    static void Exact(HNode&);
    static void Bounded(HNode&);
    static void RawBounded(HNode&);
    static void BoundedInEdge(HNode&);
    // static di BestFirst(int); 
    // static di BestFirstNoBound(); 
    static void UpdateAP(const int);
    static double MarginalAPOf(const int, const int, const double);

   

    // Lower
    static void LBLoad(char *);
//    static double LBCalculate(const vector<vector<int> >&, vector<di>&);
    static double LBCalculate(const vector<int>&);
    static int LBSelect();
    static double LBMarginalAPOf(const int, const int, const double);

    // Upper Approx
    static void UBLoad(char *);
    static double UBCalculate(double B);

    static S GreedyBound(Q, int);

    static int BarIndex(double);

public:
    static double GetKUWeight(int,int);
    static vector<IL> il_;
    static double WeightedScore(vector<int> seeds,Q q);
    static vector<double> GetItem(vector<int> );
    static void BuildKeywordNeighborhoodBound(char *);
    static void BuildNeighborhoodBound(char *);
    static void BuildInvertedList(char *);
    static void Init();
    static void Offline(char *);
    static S BestEffort(Q);
    static S Neighborhood(Q q);



    // load
    static void LoadTIC(char *,char *);
    static double theta_;
    static double theta2_;


    // query
    static void InitGreedy(char *);
    static void InitGreedyInEdge(char *);
    static S Greedy(Q);
    static S Inflex(Q);
    static S GreedyRawBound(Q);
    static S GreedyInEdge(Q);
    static S GreedyNoBound(Q);
    static void InitApprox(char *, char *, char *);
    static S Approx(Q, double epsilon);

//    static vector<vector<int> > LBGenerate(vector<int>);

    static void Reset();
    static void PrintArguments();
};

#endif
