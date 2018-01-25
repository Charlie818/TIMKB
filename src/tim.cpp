#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>    // std::max
#include <time.h>    // std::max
#include <utility>
#include <map>    // std::max

#include "tim.h"

#define INF 2147483647
#define MIN 0.0001
// #define MIN 0.0001

double* TIM::ap_;
double* TIM::tap_;
vector<HNode> TIM::H_;
Q TIM::q_;
double TIM::theta_;
double TIM::theta2_;
int TIM::round_;
double *TIM::dist_;
int *TIM::seen_;
int *TIM::seen_idx_;
int *TIM::children_;
int *TIM::pred_;
double *TIM::bars_;
vector<double*> TIM::hist_;

vector<vector<di> > TIM::iees_cache_;


vector<di> TIM::L_;
// vector<vector<di> > TIM::USet_;
vector<vector<double> > TIM::max_topic_;
// double *TIM::max_in_edge_;

int TIM::n_;
bool *TIM::used_;
double *TIM::lap_;
vector<UBSample> TIM::ubsamples_;
vector<LBSample> TIM::lbsamples_;


double exactTimer; clock_t eStart;
double boundTimer; clock_t bStart;
double dTimer; clock_t dStart; double dSize;
double apTimer; clock_t apStart;
double pTimer; clock_t pStart;
double BestFirstTimer;clock_t BStart;
double iterTimer;clock_t itStart;
int boundCounter;
int exactCounter;
int dCounter;
int insertCounter;
clock_t start;
double tTimer; clock_t tStart;
vector<IL> TIM::il_;
vector<vector<NB> > TIM::nb_;
vector<NB> TIM::nb2_;
map<int,int> cursorL;
di cursor;//pointer to the L with max benefit and indicate the sum benefit of all L
di cursor2;
double cost;

map<int,bool> uid_used;

double TIM::WeightedScore(vector<int> seeds,Q q){
    srand (time(NULL));
    int n_times = 100;

    int h, t, k_top = seeds.size();
    vector<int> list;
    set<int> active;
    double resultSize = 0;
    for (int it=0; it < n_times; it++) {
        list.clear();
        active.clear();
        for (int i = 0; i < k_top; i++) 
        {
            list.push_back(seeds[i]);
            active.insert(seeds[i]);
            resultSize+=GetKUWeights(seeds[i],q.keywords);
        }
        // resultSize += k_top;

        h = 0;
        t = k_top;

        while (h<t) {
            int k = TIC::GetOutNeighbor(list[h]);
            for (int i=0; i<k; i++) {
                Edge e = TIC::GetOutEdge(list[h], i);
                if (active.count(e.v) == 1) continue;
                if (((double)rand()/(double)RAND_MAX) < TIC::PProb(e.w, q.item)) {
                    list.push_back(e.v);
                    active.insert(e.v);
                    t++;
                    resultSize+=TIM::GetKUWeights(TIC::GetName(e.v),q.keywords);
                }
            }
            h++;
        }
    }
    return (double)resultSize / (double)n_times;
}

bool check_used(int uid){
    map<int ,bool >::iterator l_it;;
    l_it=uid_used.find(uid);
    if(l_it==uid_used.end()){
        printf("check used %d not find\n",uid );
        return 0;
    }
    else
        return uid_used[uid];
}

bool tim_L_max_cmp(const di&a, const di&b) {
    return a.first > b.first;
}

bool max_candidate_cmp(const HNode&a, const HNode&b) {
    return a.benefit < b.benefit;
}

bool sort_by_max_ub(const NB &nb1,const NB &nb2){
    return nb1.max_ub>nb2.max_ub;
}
void TIM::Reset() {
    #ifdef VERBOSE
    printf("Reseting TIM workspace ... Reset is needed for each query.\n");
    #endif

    TIC::resetCache();

    round_ = 0;

    H_.clear();
    
    cost=0.0;

   // Iee_.resize(n_);
    iees_cache_.resize(n_);

   // pre_.resize(n_);
   // In_.resize(n_);

    for (unsigned int i = 0; i < n_; i++) ap_[i] = 0;
    for (unsigned int i = 0; i < n_; i++) tap_[i] = 0;
    for(int i=0;i<n_;i++)uid_used[i]=false;
    // reset timer 
    iterTimer = 0;
    BestFirstTimer = 0;
    insertCounter = 0;
    exactCounter = 0;
    boundCounter = 0;
    dCounter = 0;
    exactTimer = 0;
    boundTimer = 0;
    dTimer = 0;
    apTimer = 0;
    pTimer = 0;
    dSize = 0;
    tTimer = 0;  // test
}

void TIM::Init() {

    L_.clear();
    n_ = TIC::GetN();

    ap_ = new double[n_];
    tap_ = new double[n_];
    used_ = new bool[n_];
    for (unsigned int i = 0; i < n_; i++) used_[i] = false;

    dist_ = new double[n_];
    seen_ = new int[n_];
    seen_idx_ = new int[n_];
    children_ = new int[n_];
    pred_ = new int[n_];
   // is_border_ = new bool[n_];
    for (unsigned int i = 0; i < n_; i++) dist_[i] = INF;
    for (unsigned int i = 0; i < n_; i++) seen_idx_[i] = n_;

    #ifdef VERBOSE
    PrintArguments();
    #endif
}

void TIM::PrintArguments() {
    printf("N:%d M:%d Z:%d theta:%lg theta2 %lg\n", n_, TIC::GetM(), TIC::GetZ(), theta_, theta2_);
}

void TIM::InitGreedy(char* filename) {
    Init();
    #ifdef VERBOSE
    printf("Initializing Greedy BestFirst framework...\n");
    #endif

    // Load \List
    FILE *f = fopen(filename, "r");
    int n, id, size; double d, d2;
    fscanf(f, "%d", &n);
    assert(n == n_);
    
    L_.resize(n);
    hist_.resize(n);
    bars_ = new double[NBARS];
    double h=log(theta_)/NBARS;
    for (int b = 0; b < NBARS; b++) bars_[b] = exp((b+1)*h);
    for (int i = 0; i < n; i++) hist_[i] = new double[NBARS];
    for (int i = 0; i < n; i++) {
        fscanf(f, "%d\t%lg\t%d", &id, &d, &size);
        L_[i] = di(d, id);
        for (int b = 0; b < NBARS; b++)
            fscanf(f, "\t%lg", &(hist_[id][b]));
    }

    fclose(f);
}

double TIM::GetIDF(int kid){
    return il_[kid].idf;
}
double TIM::GetIDF(vector<int> keywords){
    double res=0.0;
    for(int i=0;i<keywords.size();i++)
        res+=il_[keywords[i]].idf;
    return res;
}

double TIM::GetKUWeight(int uid,int kid){
    // printf("%d %d %d\n",uid,kid,TIC::GetInterest(uid,kid) );
    if(TIC::GetInterest(uid,kid)<0)return 0.0;
    vector<tuple> list=il_[kid].list;
    for(int i=0;i<list.size();i++){
        if(list[i].first==uid)return list[i].second;
    }
    return 0.0;
}
double TIM::GetKUWeights(int uid,vector<int> kids){
    double res=0.0;
    for(int i=0;i<kids.size();i++){
        // printf("uid %d kid %d res %lf \n",uid,kids[i],GetKUWeight(uid,kids[i]) );
        res+=GetKUWeight(uid,kids[i]);
    }
    return res;
}
//Bayesian Transition
vector<double> TIM::GetItem(vector<int>keywords){
    vector<double> item;
    item.resize(TIC::GetZ());
    double bottom=0.0;
    
    // for(int i =0;i<keywords.size();i++){
    //     printf("keyword %d",keywords[i]);
    //     for(int z=0;z<TIC::GetZ();z++){
    //         printf("%lf ",TIC::GetTopic(keywords[i])[z] );
    //     }
    //     printf("\n");
    // }
    // for(int i=0;i<TIC::GetZ();i++){
    //     printf("%.4lf ",item[i] );
    // }
    // printf("\n");
    for(int i=0;i<TIC::GetZ();i++){
        double tmp=1.0;
        double z=TIC::GetTopic(-1)[i];
        for(int j=0;j<keywords.size();j++){
            tmp*=TIC::GetTopic(keywords[j])[i];
        }
        bottom+=tmp*z;
    }
    // printf("bottom %lf\n",bottom );
    for(int i=0;i<TIC::GetZ();i++){
        double tmp=1.0;
        double z=TIC::GetTopic(-1)[i];
        for(int j=0;j<keywords.size();j++){
            tmp*=TIC::GetTopic(keywords[j])[i];
        }
        item[i]=tmp*z/bottom;
    }
    printf("item: ");
    for(int i=0;i<TIC::GetZ();i++){
        printf("%.4lf ",item[i] );
    }
    printf("\n");
    return item;
}

//init the cursor for query keywords and line up L and sum up the top one
void TIM::CursorInit(Q q){
    vector<int> keywords=q.keywords;
    cursorL.clear();
    cursor.first=0.0;
    cursor.second=-1;
    for (int i = 0; i < keywords.size(); ++i){
        int kid=keywords[i];
        if(nb_[kid].size()==0){
            cursorL[kid]=-1;
            continue;
        }
        cursorL[kid]=0;
        if(nb_[kid][0].max_ub>cursor.first){
            cursor.first=nb_[kid][0].max_ub;
            cursor.second=kid;
        }
    }
    cursor.first=0.0;
    for (int i = 0; i < keywords.size(); ++i)
    {
        int kid=keywords[i];
        if(nb_[kid].size()==0)continue;
        cursor.first+=nb_[kid][0].max_ub;
    }
    // printf("cursor ub %lf id %d\n",cursor.first,cursor.second );
}

void TIM::BuildInvertedList(char * filename){
    FILE *f = fopen(filename, "r");
    int keyword_num;
    fscanf(f,"%d",&keyword_num);
    printf("%d %d\n",keyword_num,TIC::GetK() );
    assert(keyword_num==TIC::GetK());
    for (int i = 0; i < keyword_num; ++i)
    {
        IL il;
        vector<tuple> l;
        int users,kid;
        double idf;
        fscanf(f,"%d %lf %d",&kid,&idf,&users);
        for (int j = 0; j < users; ++j)
        {
            int uid;
            double weight;
            fscanf(f,"%d %lf",&uid,&weight);
            l.push_back(tuple(uid,weight));
        }
        il.kid=kid;
        il.list=l;
        il.idf=idf;
        // printf("kid %d users:%d size %d\n",i,users,l.size() );
        il_.push_back(il);
    }
    fclose(f);
}

void TIM::BuildNeighborhoodBound(char * filename){
    FILE *f = fopen(filename, "r");
    int n,Z;
    fscanf(f,"%d %d",&n,&Z);
    assert(Z==TIC::GetZ());
    nb2_.resize(n);
    for (int i = 0; i < n; ++i)
    {
        int uid;
        double max_ub;
        fscanf(f,"%d %lf",&uid,&max_ub);
        vector<double> ub;
        for (int z = 0; z < Z; z++)
        {
            double weight;
            fscanf(f,"%lf",&weight);
            ub.push_back(weight);
        }
        NB t;
        t.uid=uid;
        t.ub=ub;
        t.max_ub=max_ub;
        nb2_[i]=t;
    }
    sort(nb2_.begin(),nb2_.end(),sort_by_max_ub);
    // for(int i=0;i<nb2_.size();i++){
    //     printf("id %d max_ub %lf\n",nb2_[i].uid,nb2_[i].max_ub );
    // }
    fclose(f);
}

void TIM::BuildKeywordNeighborhoodBound(char * filename){
    FILE *f = fopen(filename, "r");
    int keyword_num,Z;
    fscanf(f,"%d %d",&keyword_num,&Z);
    assert(Z==TIC::GetZ());
    assert(keyword_num==TIC::GetK());
    nb_.resize(keyword_num);
    for (int i = 0; i < keyword_num; ++i)
    {
        int kid,users;
        fscanf(f,"%d %d",&kid,&users);
        // printf("kid %d users %d \n",kid,users );
        vector<NB> nb;
        nb.resize(users);
        for (int j = 0; j < users; ++j)
        {
            int uid;
            double max_ub;
            fscanf(f,"%d %lf",&uid,&max_ub);
            vector<double> ub;
            for (int z = 0; z < Z; z++)
            {
                double weight;
                fscanf(f,"%lf",&weight);
                ub.push_back(weight);
            }
            NB t;
            t.uid=uid;
            t.ub=ub;
            t.max_ub=max_ub;
            // vector<double>::iterator it=max_element(begin(ub),end(ub));
            // t.max_ub=(TIM::GetKUWeight(uid,kid)+*it)/TIC::GetCost(uid);
            nb[j]=t;
        }
        sort(nb.begin(),nb.end(),sort_by_max_ub);
        nb_[i]=nb;

    }
    fclose(f);
}

void TIM::InitGreedyInEdge(char* filename) {
    Init();
    #ifdef VERBOSE
    printf("Initializing Greedy BestFirst In Edge framework...\n");
    #endif

    // Load \Uset
    vector<di> rst;
    int id, n; double influ;
    FILE *f = fopen(filename, "r");
    fscanf(f, "%d", &n);
    assert (n == n_);
    L_.resize(n);
    max_topic_.resize(n);
    for	(int i = 0; i < n_; i++) {
        fscanf(f, "%d %lg", &id, &influ);
  
        L_[i] = di(influ/TIC::GetCost(id), id);
        max_topic_[i].resize(TIC::GetZ());
        for (int z = 0; z < TIC::GetZ(); z++)
            fscanf(f, "%lg", &(max_topic_[i][z]));
    }
    sort(L_.begin(), L_.end(), tim_L_max_cmp);
    fclose(f);
}

// void TIM::InitApprox(char *gfile, char *ufile, char *lfile) {
//     Offline(gfile);
//     //    InitGreedy(gfile);
//     #ifdef VERBOSE
//     printf("Initializing Approximation framework...\n");
//     #endif
//     UBLoad(ufile);
//     LBLoad(lfile);
//     #ifdef VERBOSE
//     printf("Initializing Approximation framework done.\n");
//     #endif
// }
//Besteffort
S TIM::BestEffort(Q q){
    start = clock();
    Reset();
    CursorInit(q);
    q_ = q;
    di s;
    S S2;
    while(true){
        BStart=clock();
        s = BestFirstWithCost(q);
        BestFirstTimer+=double(clock() - BStart) / CLOCKS_PER_SEC;
        if(s.second==-1)break;
        S2.ids.push_back(s.second);
        S2.minflu.push_back(s.first);
        S2.cost.push_back(TIC::GetCost(s.second));
        // printf("select id:%d influ:%lf cost:%lf ratio:%lf\n",s.second,s.first,TIC::GetCost(s.second),s.first/TIC::GetCost(s.second));
        UpdateAP(TIC::GetId(s.second));
        used_[TIC::GetId(s.second)] = true; 
    }
    for (int i = 0; i < S2.ids.size(); i++) {
        used_[TIC::GetId(S2.ids[i])] = false;
    }
    // double S2_minflu=0.0;
    // for(int i=0;i<S2.ids.size();i++){
    //     int id=S2.ids[i];
    //     printf("name: %d cost %lf benefit %lf\n",id,TIC::GetCost(id),S2.minflu[i]);
    // }
    // #ifdef VERBOSE
    // printf("total cost %lf\n",cost );
    // printf("influ:S1:%lf S2:%lf\n",S1.minflu[0],S2_minflu );
    // printf("bound timer:%lf exact timer:%lf bound cnt: %d exact cnt: %d \n",boundTimer/CLOCKS_PER_SEC,exactTimer/CLOCKS_PER_SEC ,boundCounter,exactCounter );
    // printf("S2 consuming time: %lg \n", double(clock() - start) / CLOCKS_PER_SEC);
    // printf("BestFirstTimer:%lg\n",BestFirstTimer );
    // printf("insert cnt %d iterTimer %lg \n",insertCounter,iterTimer );
    // #endif
    return S2; 
}
void TIM::UpdateAP(const int u) {
    int v;
    vector<di>& rst = iees_cache_[u];
   // printf("round %d, (%d, %lg), size %d\n", round_, u, theta_, rst.size());
    for (unsigned int i = 0; i < rst.size(); i++) {
        v = rst[i].second;
        ap_[v] += rst[i].first;
        tap_[v] = ap_[v];
       // In_[v].insert(rst[pre[i]].second);
    }
}

double TIM::MarginalAPOf(const int v, const int p, const double w) {
   // apStart = clock();

    if (used_[v]) return 0;
   // double ap_v = 1.0;
   // for (unsigned int i = 0; i < TIC::GetInNeighbor(v); i++) {
   //     Edge e = TIC::GetInEdge(v, i);
   //     if (tap_[e.u] > 0) {
   //         ap_v *= (1 - TIC::GetInEdgeWeight(v, i, q_.item)*tap_[e.u]);
   //     }
   // }
   // printf("In_[%d]={", v);
   // for (set<int>::iterator it = In_[v].begin(); it != In_[v].end(); it++) {
   //     printf("%d ", *it);
   // }
   // printf("}\n");
   // ap_v = 1.0 - ap_v;
    double ap_v = ap_[v], apw = ap_[p]*w == 1 ? 0.99 : ap_[p]*w;
    ap_v = 1.0 - ap_v;
    ap_v /= (1.0 - apw);
    ap_v *= (1.0 - tap_[p]*w);
    ap_v = 1.0 - ap_v;
    tap_[v] = ap_v;
    #ifdef VERBOSE
    if (ap_v < ap_[v]) {
        printf("ap(%d)=%lg, tap=%lg\n", v, ap_[v], tap_[v]);
    }
   // assert (ap_v >= ap_[v]);
    #endif
    
   // apTimer += double(clock() - apStart);
    return ap_v - ap_[v];
}

int TIM::Dijkstra(int u, double max) {
    dist_[u] = 0;
    int top = 0;
    int bottom = 0;
    seen_[top++] = u;
    children_[bottom++] = u;
    seen_idx_[u] = -1;

   // rst.push_back(di(1.0f, u));
   // pre.push_back(0);

    // variables used in the while loop.
    double w_dist; int v, w;
    double tmp;

    while (top > 0) {
        v = seen_[0];
       // printf("Dijkstra pop (%d, %lg) %lf\n", v, exp(-dist_[v]),exp(-max));

        if (dist_[v] < max) {
            seen_idx_[v] = -1;
        } else break;

        // shortest distance vertex for sure.
        if (v != u) {
           // rst.push_back(di(exp(-dist_[v]), v));
            children_[bottom++] = v;
           // is_border_[pred_[v]] = false;
        }

        for (unsigned int i = 0; i < TIC::GetOutNeighbor(v); i++) {
            Edge e = TIC::GetOutEdge(v, i);
            // printf("e(%d, %d) %lf %lf\n", e.u, e.v,e.w[0],e.w[1]);
            w = e.v;

            if (used_[w] || seen_idx_[w] < 0) continue;

           // pStart = clock();
            w_dist = -log(TIC::GetOutEdgeWeight(v, i, q_.item));
           // pTimer += double(clock() - pStart);

           // tStart = clock();
            if (w_dist + dist_[v] < dist_[w]) {
                dist_[w] = w_dist + dist_[v];
                pred_[w] = v;
                // push
                int j;
                if (seen_idx_[w] >= n_) {  // not seen before
                    seen_[top] = w;
                    j = top++;
                } else {
                    j = seen_idx_[w];
                }
                int x = (j-1)/2;
                tmp = dist_[seen_[j]];
                while (j > 0) {
                    if (dist_[seen_[x]] > tmp) {
                        seen_[j] = seen_[x];
                        if (seen_idx_[seen_[j]] < n_) seen_idx_[seen_[j]] = j;
                        j = x; x = (j-1)/2;
                    } else break;
                }
                seen_[j] = w;
                seen_idx_[w] = j;
            }
           // tTimer += double(clock() - tStart);
        }

        // pop
        seen_[0] = seen_[--top];
        if (!top) break;
        int j = 0, x = j*2 + 1;
        tmp = dist_[seen_[j]];
        while (x < top) {
            if (x+1 < top && dist_[seen_[x+1]] < dist_[seen_[x]]) x=x+1;
            if (dist_[seen_[x]] < tmp) {
                seen_[j] = seen_[x];
                seen_idx_[seen_[j]] = j;
                j = x; x = j*2+1;
            } else break;
        }
        seen_[j] = seen_[top];
        if (seen_idx_[seen_[j]] < n_) seen_idx_[seen_[j]] = j;
    }

    // restore
    for (unsigned int i = 0; i < bottom; i++) {
       // dist_[rst[i].second] = INF;
       // seen_idx_[rst[i].second] = n_;
       // dist_[children_[i]] = INF;
        seen_idx_[children_[i]] = n_;
    }

    for (unsigned int i = 0; i < top; i++) {
        dist_[seen_[i]] = INF;
        seen_idx_[seen_[i]] = n_;
    }

   // dTimer += double(clock() - dStart);

   // return rst;
   // dSize += bottom;
    return bottom;
}

HNode TIM::pop_H(){
    HNode t=H_.front();
    pop_heap(H_.begin(),H_.end(), max_candidate_cmp); 
    H_.pop_back();
    return t;
}

void TIM::push_H(HNode t){
    H_.push_back(t); 
    push_heap(H_.begin(), H_.end(), max_candidate_cmp);
}

di TIM::BestFirstWithCost(Q q) {
    HNode u;
    do{
        itStart = clock();
        InsertCandidatesWithCost(q);
        if(H_.empty())return di(-1,-1);
        u=pop_H();
        // printf("u id %d cost %lf \n",u.name,u.cost );
        if(u.cost+cost>q.B){
            iterTimer+=double(clock() - itStart) / CLOCKS_PER_SEC;
            continue;
        }
        // printf("top u:%d cost:%lf benefit %lf status %d  H:size %d\n",u.name,u.cost,u.benefit,u.status,H_.size()+1);
        if(u.status==INITIAL){
            u.benefit=EstMarginUB(u,q)/u.cost;
            u.status=ESTIMATED;
            push_H(u);
        }else if(u.status==ESTIMATED){
            if(u.round==round_){
                double exact=CalcMargin(u,q);
                if(u.benefit<exact/u.cost){
                    printf("error id %d %lf %lf\n",u.id,u.benefit,exact/u.cost);
                }
                u.influ=exact;
                u.benefit=exact/u.cost;
                u.status=COMPUTED;
                push_H(u);
            }
            else{
                u.influ=EstMarginUB(u,q);
                u.benefit=u.influ/u.cost;
                u.status=ESTIMATED;
                u.round=round_;
                push_H(u);
            }
        }else if(u.status==COMPUTED){
            if(u.round==round_){
                cost+=u.cost;
                round_++;
                printf("select name %d influ %lf cost%lf benefit%lf\n",u.name,u.influ,u.cost,u.benefit );
                return di(u.benefit*u.cost,u.name);
            }else{
                u.influ=CalcMargin(u,q);
                u.benefit=u.influ/u.cost;
                u.status=COMPUTED;
                u.round=round_;
                push_H(u);
            }
        }
    }while(q.B-cost>=1);
    return di(-1,-1);
}

void TIM::InsertCandidatesWithCost(Q q){
    int kid=cursor.second;
    double sum=cursor.first;
    double maxH = H_.empty() ? 0 : H_.front().benefit;
    while(sum>maxH){
        if(kid==-1)break;
        // printf("kid %d sum %lf idx %d uid %d \n",kid,sum,cursorL[kid],nb_[kid][cursorL[kid]].uid );
        vector<NB> nbs = nb_[kid];
        int idx=cursorL[kid];
        int uid=nbs[idx].uid;
        double maxL=nbs[idx].max_ub;
        // printf("kid %d uid %d\n",kid,uid );
        HNode c;
        c.id =  TIC::GetId(uid);
        c.name = uid;
        c.status = ESTIMATED;
        c.cost=TIC::GetCost(c.name);
        c.benefit = EstMarginUB(c,q)/c.cost;
        c.round=round_;
        insertCounter++;
        push_H(c);
        uid_used[c.name]=true;
        // printf("insert name %d cost %lf benefit %lf H_ size %d\n",c.name,c.cost,c.benefit,H_.size() );

        // //TODO
        // //update sum
        // sum-=maxL;
        // int offset=1;
        // int next_id=nb_[kid][idx+offset].uid; 
        // while((offset+idx)<nb_[kid].size()&&check_used(next_id)){
        //     offset++;
        //     next_id=nb_[kid][idx+offset].uid; 
        // }
        // if((offset+idx)>=nb_[kid].size())cursorL[kid]=-1;
        // else{
        //     cursorL[kid]=idx+offset;
        //     sum+=nbs[idx+offset].max_ub;
        // }
                

        double value=0.0;
        for(map<int,int>::iterator i =cursorL.begin();i!=cursorL.end();i++){
            if(i->second==-1)continue;
            int kid=i->first,idx=i->second;
            while(idx<nb_[kid].size()&&check_used(nb_[kid][idx].uid))
                idx++;
            
            if(idx==nb_[kid].size()){
                cursorL[kid]=-1;
                continue;
            }
            cursorL[kid]=idx;
            // printf("check %d %d %d\n",idx_of_L,idx_in_L,check_used(nb_[idx_of_L][idx_in_L].uid));
            if(nb_[kid][idx].max_ub>value){
                cursor.second=kid;
                value=nb_[kid][idx].max_ub;
            }
        }
        cursor.first=0.0;
        for (int i = 0; i < q.keywords.size(); ++i)
        {
            int kid=q.keywords[i];
            if(cursorL[kid]==-1)continue;
            cursor.first+=nb_[kid][cursorL[kid]].max_ub;
        }
        kid=cursor.second;
        sum=cursor.first;
        maxH = H_.front().benefit;
        // printf("cursor 0 %d cursor 1 %d\n",cursorL[0],cursorL[1] );
        // printf("sum %lf maxH %lf cursor %d\n",sum,maxH,kid);
    }
}

double TIM::EstMarginUB(HNode seed,Q q){
    if(seed.id>TIC::GetN())printf("error\n");
    bStart = clock();
    boundCounter++;
    
    //calc the center node benefits wrt keywords
    double inflUB=GetKUWeights(seed.name,q.keywords);
    // printf("name %d influUB %lf\n",seed.name,inflUB);

    for (int i = 0; i < q.keywords.size(); ++i)
    {
        int kid =q.keywords[i];
        for (int j=cursorL[kid];j<nb_[kid].size();j++){
            NB t=nb_[kid][j];
            if(t.uid!=seed.name)continue;
            for (int z = 0; z < TIC::GetZ(); ++z)
            {
                inflUB+=t.ub[z]*q.item[z];
            }
        } 
    }
    // printf("finish\n");
    inflUB *= 1 - ap_[seed.id];
    // printf("inflUB:%lf\n",inflUB );
    boundTimer += double(clock() - bStart);
    return inflUB;
}

double TIM::CalcMargin(HNode seed,Q q){
    eStart = clock();
    exactCounter++;

    //?
    int bottom = Dijkstra(seed.id, -log(theta_));
    // printf("bottom %d \n",bottom );
    /* calculate marginal influence of seed */
    tap_[seed.id] =1.0;
    double rt = (tap_[seed.id] - ap_[seed.id])*GetKUWeights(seed.name,q.keywords), minflu;
    // printf("exact id %d init %lf \n",seed.name,rt/TIC::GetCost(seed.name) );
    iees_cache_[seed.id].resize(bottom);
    iees_cache_[seed.id][0] = di(tap_[seed.id] - ap_[seed.id], seed.id);

    int id, p;
    for (unsigned int i = 1; i < bottom; i++) {
        id = children_[i]; p = pred_[id];
        
        minflu = MarginalAPOf(id, p, exp(dist_[p]-dist_[id]));
        rt += minflu*GetKUWeights(TIC::GetName(id),q.keywords);
        iees_cache_[seed.id][i] = di(minflu, id);
    }
    
    // printf("exact id %d marginal %lf \n",seed.name,rt/TIC::GetCost(seed.name) );

    for (int i = 0; i < bottom; i++) {
        /* restore */
        dist_[children_[i]] = INF;
        tap_[children_[i]] = ap_[children_[i]];
    }
    exactTimer += double(clock() - eStart);
    return rt;
}

S TIM::Neighborhood(Q q){
    start = clock();
    Reset();
    cursor2=di(nb2_[0].max_ub,0);
    q_ = q;
    di s;
    S S2;
    while(true){
        BStart=clock();
        s = NeighborhoodWithCost(q);
        BestFirstTimer+=double(clock() - BStart) / CLOCKS_PER_SEC;
        if(s.second==-1)break;
        S2.ids.push_back(s.second);
        S2.minflu.push_back(s.first);
        S2.cost.push_back(TIC::GetCost(s.second));
        // printf("select id:%d influ:%lf cost:%lf ratio:%lf\n",s.second,s.first,TIC::GetCost(s.second),s.first/TIC::GetCost(s.second));
        UpdateAP(TIC::GetId(s.second));
        used_[TIC::GetId(s.second)] = true; 
    }
    for (int i = 0; i < S2.ids.size(); i++) {
        used_[TIC::GetId(S2.ids[i])] = false;
    }
    // double S2_minflu=0.0;
    for(int i=0;i<S2.ids.size();i++){
        int id=S2.ids[i];
        printf("name: %d cost %lf benefit %lf\n",id,TIC::GetCost(id),S2.minflu[i]);
    }
    return S2; 
}

di TIM::NeighborhoodWithCost(Q q) {
    HNode u;
    do{
        itStart = clock();
        InsertCandidates2(q);
        if(H_.empty())return di(-1,-1);
        u=pop_H();
        // printf("u id %d cost %lf \n",u.name,u.cost );
        if(u.cost+cost>q.B){
            iterTimer+=double(clock() - itStart) / CLOCKS_PER_SEC;
            continue;
        }
        // printf("top u:%d cost:%lf benefit %lf status %d  H:size %d\n",u.name,u.cost,u.benefit,u.status,H_.size()+1);
        if(u.status==ESTIMATED){
            if(u.round==round_){
                double exact=CalcMargin(u,q);
                if(u.benefit<exact/u.cost){
                    printf("error\n");
                }
                u.benefit=exact/u.cost;
                u.status=COMPUTED;
                push_H(u);
            }
            else{
                u.benefit=EstMarginUB2(u,q)/u.cost;
                u.status=ESTIMATED;
                u.round=round_;
                push_H(u);
            }
        }else if(u.status==COMPUTED){
            if(u.round==round_){
                cost+=u.cost;
                round_++;
                // printf("select name %d ratio %lf cost%lf \n",u.name,u.benefit,u.cost );
                return di(u.benefit*u.cost,u.name);
            }else{
                u.benefit=CalcMargin(u,q)/u.cost;
                u.status=COMPUTED;
                u.round=round_;
                push_H(u);
            }
        }
    }while(q.B-cost>=1);
    return di(-1,-1);
}

void TIM::InsertCandidates2(Q q){
    int pt=cursor2.second;
    double sum=cursor2.first;
    double maxH = H_.empty() ? 0 : H_.front().benefit;
    while(sum>maxH){
        NB nb = nb2_[pt];
        int uid=nb.uid;
        double maxL=nb.max_ub;

        HNode c;
        c.id =  TIC::GetId(uid);
        c.name = uid;
        c.status = ESTIMATED;
        c.cost=TIC::GetCost(c.name);
        c.vec=nb.ub;
        c.benefit = EstMarginUB2(c,q)/c.cost;
        c.round=round_;
        
        insertCounter++;
        push_H(c);
        // printf("insert name %d cost %lf benefit %lf H_ size %d\n",c.name,c.cost,c.benefit,H_.size() );

        pt++;
        if(pt==nb2_.size())sum=-1;
        else sum=nb2_[pt].max_ub;
        cursor2=di(sum,pt);
        maxH = H_.front().benefit;
    }
}

double TIM::EstMarginUB2(HNode seed,Q q){
    if(seed.id>TIC::GetN())printf("error\n");
    bStart = clock();
    boundCounter++;
    
    //calc the center node benefits wrt keywords
    double inflUB=GetKUWeights(seed.name,q.keywords);
    double temp=0.0;
    // printf("name %d influUB %lf\n",seed.name,inflUB);
    for (int z = 0; z < TIC::GetZ(); ++z)
    {
        temp+=seed.vec[z]*q.item[z];
    }
    inflUB+=temp*GetIDF(q.keywords);
    inflUB *= 1 - ap_[seed.id];
    boundTimer += double(clock() - bStart);
    return inflUB;
}


// void TIM::BoundedInEdge(HNode& seed) {
//     boundCounter += 1;
//     bStart = clock();
 
//     double minflu = 1;
//     for (unsigned int z = 0; z < TIC::GetZ(); z++) {
//         minflu += q_.item[z] * max_topic_[seed.id][z];
//     }
//     minflu *= 1 - ap_[seed.id];

//     #ifdef VERBOSE
//     //    printf("Bounded (%d, %lg)\n", seed.id, minflu);
//     #endif
//     seed.infl = minflu;  /* lazy forward */
//     seed.round = round_;
//     seed.status = ESTIMATED;

//     boundTimer += double(clock() - bStart);
//     #ifdef EXAMPLE
//     printf("Bounded (%d, %lg)\n", seed.id, minflu);
//     #endif
// }


//Approximation Mathod
S TIM::Approx(Q q, double epsilon) {

    Reset();
    q_ = q;

    double UB = UBCalculate(q.B);
    double LB=0.0;
    S rsp;
    di s;
    #ifdef VERBOSE
    printf("query\n");
    PrintQ(q.item);
    #endif
    lap_ = new double[n_];
    for (int i = 0; i < n_; i++) lap_[i] = 0;

    LBSample lb = lbsamples_[LBSelect()];
    #ifdef VERBOSE
    printf("lower sample\n");
    PrintQ(lb.item);
    #endif
    vector<int> candidates; //convert to idx
    double candidates_cost=0.0;
    for(int i=0;i<lb.seeds.size();i++){
        if(lb.cost[i]<=q.B)
            candidates.push_back(TIC::GetId(lb.seeds[i]));
        else{
            LB=lb.influ[i];
            candidates_cost=lb.cost[i];
            break;
        }
    }
    #ifdef VERBOSE
    printf("lb %d candidates %d UB %lf LB %lf \n",lb.seeds.size(),candidates.size(),UB,LB );
    #endif
    double   Real = 0,Real_cost=0.0;
    int select_idx = -1;  /* the index of candidate that are selected by Greedy*/
    int step = 5;
    tTimer = 0;
    if (LB > UB*epsilon) {
        for (unsigned int i = 0; i < candidates.size(); i++) {
            rsp.ids.push_back(candidates[i]);
            rsp.minflu.push_back(0);
        }

    #ifdef VERBOSE
        printf("%d, %d,  %lg, %lg, ", boundCounter, exactCounter,  UB, LB);
    #endif
        /* restore */
        for (int i = 0; i < rsp.ids.size(); i++) used_[rsp.ids[i]] = false;

        return rsp;
    }
    #ifdef VERBOSE
    printf("spread %lf\n",TIC::Spread(candidates,q.item) );
    #endif
    while(true){
        if (round_%step == 0) {
           // tStart = clock();

           LB = Real + LBCalculate(candidates);
          
           // tTimer += double(clock() - tStart);
            #ifdef VERBOSE
                printf("real %lf lb %lf \n",Real,LBCalculate(candidates) );
                printf("round %d, Ubound %lg, Lbound %lg\n", round_, UB, LB);
            #endif
        }
        if (LB > UB*epsilon) {
            for (unsigned int i = 0; i < candidates.size(); i++) {
                rsp.ids.push_back(candidates[i]);
                rsp.minflu.push_back(0);
            }

            #ifdef VERBOSE
                printf("%d, %d, %d, %lg, %lg, ", boundCounter, exactCounter, round_, UB, LB);
            
                printf("size:%d spread %lf\n",rsp.ids.size(),TIC::Spread(rsp.ids,q.item) );
            #endif
            /* restore */
            for (int i = 0; i < rsp.ids.size(); i++) used_[rsp.ids[i]] = false;
            return rsp;
        }
        s = BestFirstWithCost(q);
        if(s.second==-1)break;
        rsp.ids.push_back(TIC::GetId(s.second));
        rsp.minflu.push_back(s.first);
        Real += s.first;
        // printf("s.second %d\n",s.second );
        Real_cost+= TIC::GetCost(s.second);
        UpdateAP(TIC::GetId(s.second));
        used_[TIC::GetId(s.second)] = true;
        for (unsigned int i = 0; i < candidates.size(); i++)
            if (candidates[i] == TIC::GetId(s.second)) select_idx = i;
        if (!(select_idx < 0)) {
            for (unsigned int i = select_idx; i < candidates.size()-1; i++) {
                candidates[i] = candidates[i+1];
            }
        }
        while(Real_cost+candidates_cost>q.B){
            if(candidates.size()==0){
                #ifdef VERBOSE
                printf("cannot reach upper bound,return the best effort result for sure\n");
                #endif
                for (int i = 0; i < rsp.ids.size(); i++) used_[rsp.ids[i]] = false;
                return rsp;
                break;
            }
            int pop=candidates[candidates.size()-1];
            // printf("pop %d size%d\n",pop,candidates.size()-1);
            candidates_cost-=TIC::GetCost(TIC::GetName(pop));
            candidates.pop_back();
        }
        select_idx = -1;  /* restore */
    }


    for (int i = 0; i < rsp.ids.size(); i++) {
        // restore
        used_[rsp.ids[i]] = false;
    }

    #ifdef VERBOSE
     printf("%d, %d, %d, %lg, %lg, ", boundCounter, exactCounter, round_, UB, LB);
    #endif

    return rsp;
}
void PrintQ(vector<double> item){
    for (int i = 0; i < TIC::GetZ(); ++i)
    {
        printf("%lg\t",item[i] );
    }
    printf("\n");
}
void TIM::UBLoad(char *filename) {
    #ifdef VERBOSE
    printf("Loading Upper Bound samples ...\n");
    #endif
    int n_samples, Z;
    double B;
    FILE *f = fopen(filename, "r");
    fscanf(f, "%d\t%d\t%lf", &n_samples, &Z, &B);
    // printf("Z:%d,%d\n",Z,TIC::GetZ() );
    assert(Z == TIC::GetZ());
    ubsamples_.resize(n_samples);

    for (int i = 0; i < n_samples; i++) {
        double min_ratio=0,max_ratio,b,spread;
        vector<BNode> nodes;
        while(min_ratio!=-1){
            fscanf(f, "%lg %lg %lg %lg", &min_ratio,&max_ratio,&b,&spread);
            BNode bn;
            bn.min_ratio=min_ratio;
            bn.max_ratio=max_ratio;
            bn.B=b;
            bn.spread=spread;
            nodes.push_back(bn);
        }
        ubsamples_[i].nodes=nodes;

        ubsamples_[i].item.resize(Z);
        for (int z = 0; z < Z; z++) {
            fscanf(f, "%lg", &(ubsamples_[i].item[z]));
        }
       
    }
    fclose(f);
}

void TIM::LBLoad(char *filename) {
    #ifdef VERBOSE
    printf("Loading Lower Bound samples ...\n");
    #endif
    int n_samples, Z;
    double B;
    FILE *f = fopen(filename, "r");
    fscanf(f, "%d %d %lf", &n_samples, &Z, &B);
    assert(Z == TIC::GetZ());
    lbsamples_.resize(n_samples);

    for (int i = 0; i < n_samples; i++) {
        double total=0.0,spread=0.0;
        lbsamples_[i].item.resize(Z);
        for (int z = 0; z < Z; z++) {
            fscanf(f, "%lg", &(lbsamples_[i].item[z]));
        }
        int k=0;
        fscanf(f, "%d ", &k);
        
        for (int j = 0; j < k; ++j){
            int idx;
            double influ;
            fscanf(f, "%d %lf ", &idx,&influ);
            
                // printf("i %d idx %d influ %lf\n",i,idx,influ );
            double cost=TIC::GetCost(idx);
            total+=cost;
            
            lbsamples_[i].seeds.push_back(idx);
            lbsamples_[i].cost.push_back(total);
            lbsamples_[i].influ.push_back(influ);
        }
        
       // lbsamples_[i].lbouts = lbouts;
    }
    fclose(f);
}

double TIM::LBMarginalAPOf(const int v, const int p, const double w) {
    if (used_[v]) return 0;

    double rt;
    double ap_v = 1.0;
    for (unsigned int i = 0; i < TIC::GetInNeighbor(v); i++) {
        Edge e = TIC::GetInEdge(v, i);
        if (lap_[e.u] > 0) {
           // printf("e(%d, %d, %lg), ap(%d)=%lg\n", e.u, e.v, TIC::PProb(q_.item, e.w), e.u, tap_[e.u]);
            ap_v *= (1 - TIC::GetInEdgeWeight(v, i, q_.item)*lap_[e.u]);
        }
    }
    ap_v = 1.0 - ap_v;
    rt = ap_v - lap_[v];
    lap_[v] = ap_v;

   // assert (ap_v >= ap_[v]);
   // printf("ap(%d)=%lg, tap=%lg\n", v, ap_[v], tap_[v]);
    return rt;
}


double TIM::LBCalculate(const vector<int>& seeds) {
    double rt = 0, m_ap, minflu;
    for (unsigned int i = 0; i < seeds.size(); i++) {
        minflu = 1 - lap_[seeds[i]];
        lap_[seeds[i]] = 1.0f;  // seed
        used_[seeds[i]] = true;

        int bottom = Dijkstra(seeds[i], -log(theta2_));
        // printf("bottom %d theta %lf\n",bottom,theta2_);
        int id, p;
        for (int i = 1; i < bottom; i++) {
            id = children_[i]; p = pred_[id];
            m_ap = LBMarginalAPOf(id, p, exp(dist_[p]-dist_[id]));
            minflu += m_ap;
        }
        rt += minflu;

    // #ifdef VERBOSE
    //     printf("%d lb minflu %lg\n", seeds[i], minflu);
    // #endif

        /* restore */
        for (int i = 0; i < bottom; i++) dist_[children_[i]] = INF;
    }

    for (unsigned int i = 0; i < seeds.size(); i++) {
        /* restore */
        used_[seeds[i]] = false;
    }

    for (int i = 0; i < n_; i++) lap_[i] = ap_[i];

    return rt;
}
//v1 is q   Question
double Cosin(const vector<double>& v1, const vector<double>& v2) {
    double inner_prod = 0.0f;
    int n = v1.size();
    for (int i = 0; i < n; i++) {
        if(v2[i]-v1[i]>MIN)return 0.0;
        inner_prod += v1[i]*v2[i];
    }

    double v1_norm2 = 0.0f;
    double v2_norm2 = 0.0f;
    for (int i = 0; i < n; i++) {
        v1_norm2 += v1[i]*v1[i];
        v2_norm2 += v2[i]*v2[i];
    }
    v1_norm2 = sqrt(v1_norm2);
    v2_norm2 = sqrt(v2_norm2);

    return inner_prod / (v1_norm2*v2_norm2);
}

double ParetolDistance(const vector<double>& Q, const vector<double>& U) {
    int n = Q.size();
    double rt = 0;

    for (int i = 0; i < n; i++) {
        if (Q[i] - U[i] > MIN) {
            rt = INF; break;
        }
        rt += fabs(U[i] - Q[i]);
    }

    return rt;
}


double TIM::UBCalculate(double B) {
    #ifdef VERBOSE
    printf("ub calc\n");
    #endif
    double rt;

    // KLD
    double m = TIC::GetZ(), val; int idx = 0;
    // printf("size %d\n",ubsamples_.size() );
    for (unsigned int i = 0; i < ubsamples_.size(); i++) {
        UBSample& ub = ubsamples_[i];
        // for(int j=0;j<TIC::GetZ();j++)
        //     printf("q %lf ub %lf\n",q_.item[j],ub.item[j] );
        val = ParetolDistance(q_.item, ub.item);
        // printf("KLD %lg\n", val);
        if (val < m) {
            m = val;
            idx = i;
        }
    }
    #ifdef VERBOSE
    printf("idx %d m %lf\n",idx,m);
    PrintQ(ubsamples_[idx].item);
    #endif
    int i,j=0;
    double min_ratio=0,max_ratio=0,B1=0,B2=0,spread1=0,spread2=0;
    while(B>ubsamples_[idx].nodes[j].B){
        j++;
        if(j==ubsamples_[idx].nodes.size()-1){
            i=j;
            j=-1;
            break;
        }
    }
    if(B==ubsamples_[idx].nodes[j].B){
        rt=ubsamples_[idx].nodes[j].spread;
    }else if(j==-1){
        max_ratio=ubsamples_[idx].nodes[i].max_ratio;
        B1=ubsamples_[idx].nodes[i].B;
        spread1=ubsamples_[idx].nodes[i].spread;
        rt=spread1+(B-B1)*max_ratio;
    }else if(i==-1){
        min_ratio=ubsamples_[idx].nodes[j].min_ratio;
        B2=ubsamples_[idx].nodes[j].B;
        spread2=ubsamples_[idx].nodes[j].spread;
        rt=spread2-(B2-B)*min_ratio;
        max_ratio=-1;
    }else{
        i=j-1;
        max_ratio=ubsamples_[idx].nodes[i].max_ratio;
        min_ratio=ubsamples_[idx].nodes[j].min_ratio;
        B1=ubsamples_[idx].nodes[i].B;
        spread1=ubsamples_[idx].nodes[i].spread;
        B2=ubsamples_[idx].nodes[j].B;
        spread2=ubsamples_[idx].nodes[j].spread;
        rt=min(spread1+(B-B1)*max_ratio,spread2-(B2-B)*min_ratio);
    }

    // printf("min_ratio %lf max_ratio%lf B %lf B1 %lf B2 %lf ",min_ratio,max_ratio,B,B1,B2 );
    // printf("spread1 %lf spread2 %lf rt %lf\n",spread1,spread2,rt );

    // rt=CalcUpperBound(idx,B);
    // #ifdef VERBOSE
    // printf("chosen %dth upper bound item sample:\n", idx);
    // for (unsigned int i = 0; i < ubsamples_[idx].item.size(); i++) {
    //     printf("%lg ", ubsamples_[idx].item[i]);
    // }
    // printf("bound: %lg", rt);
    // printf("\n");
    // #endif

    // Cosin
    return rt;
}

int TIM::LBSelect() {
    // KLD
    double m = 0.0f, val; int key = 0;
    for (unsigned int i = 0; i < lbsamples_.size(); i++) {
        LBSample& lb = lbsamples_[i];
        // printf("lb item %d", lb.item.size());
        val = Cosin(q_.item, lb.item);
        // printf("val %lf\n",val );
        if (val > m) {
            m = val;
            key = i;
        }
    }
    #ifdef VERBOSE
    if(m==0.0)printf("still can not find a lower sample\n");
    printf("key %d m %lf\n",key,m );
    #endif
    return key;
}

int TIM::BarIndex(double t) {
    for (int b = 0; b < NBARS; b++) {
	if (bars_[b] < t)
            return b;
    }
    return NBARS-1;
}


// //Greedy Method
// S TIM::Greedy(Q q){
//     start = clock();
//     Reset();
//     #ifdef VERBOSE
//         printf("reset takes time %lg s.\n", double(clock() - start) / CLOCKS_PER_SEC);
//     #endif
//     q_ = q;
//     di s;
//     InsertCandidatesWithCostGreedy(L2);
//     S S;
//     while(true){

//         BStart=clock();
//         s = GreedyWithCost(q);
//         BestFirstTimer+=double(clock() - BStart) / CLOCKS_PER_SEC;
//         if(s.second==-1)break;
//         S.ids.push_back(s.second);
//         S.minflu.push_back(s.first);
//         S.cost.push_back(TIC::GetCost(s.second));
//         UpdateAP(TIC::GetId(s.second));
//         used_[TIC::GetId(s.second)] = true; 
//     }
//     for (int i = 0; i < S.ids.size(); i++) {
//         used_[TIC::GetId(S.ids[i])] = false;
//     }
//     return S;
// }
// di TIM::GreedyWithCost(Q q) {
//     HNode u;
//     itStart = clock();
//     int key_idx=-1;
//     double key_infl=0.0,key_cost=0.0;
//     for(int i=0;i<H_.size();i++){
        
//         if(q.B-cost<1)continue;
//         if(used_[i])continue;
//         u=H_[i];

//         if(u.cost+cost>q.B)continue;
//         u.infl=CalcMargin(u,q)/u.cost;
//         if(u.infl>key_infl){
//             key_infl=u.infl;
//             key_idx=u.name;
//             key_cost=u.cost;
//         }
//     }
    
//     if(key_idx!=-1)cost+=key_cost;
//     // printf("key id %d key infl %lf total cost %lf\n",key_idx,key_infl,cost );
//     return di(key_infl,key_idx);
// }


