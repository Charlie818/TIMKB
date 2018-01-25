// offline methods
#include <string>
#include <stdio.h>
#include <math.h>
#include <algorithm>    // std::max
#include <sys/time.h>    // std::max
#include <utility>
#include <map>    // std::max

#include "tim.h"

bool max_lnode_cmp(const LNode&a, const LNode&b) {
    return a.infl > b.infl;
}

void USET(double theta, char *filename) {
    FILE *f = fopen(filename, "w");
    fprintf(f, "%d %d\n", TIC::GetN(), TIC::GetZ());
    double *max_u_z = new double[TIC::GetZ()];
    double *max_v_z = new double[TIC::GetZ()];
    double *influ = new double[TIC::GetZ()];
    for (int i = 0; i < TIC::GetN(); i++) {
        if (i % 10000 == 0) printf("processing vertex %d...\n", i);
        for (int z = 0; z < TIC::GetZ(); z++) influ[z] = 0;
        max_u_z = TIC::GetOutEdgeBound(i);

        vector<di> rst = TIC::Dijkstra(i, -log(theta));
        double rt = 1;
        for (unsigned int i = 1; i < rst.size(); i++) {
            int v = rst[i].second;
            int flag = 0; //if 1,then v is near u
            for (int z = 0; z < TIC::GetZ(); z++) max_v_z[z] = 0;
            for (int x = 0; x < TIC::GetInNeighbor(v); x++) {
                Edge e = TIC::GetInEdge(v, x);
                if (i == e.u)flag = 1;
                for (int z = 0; z < TIC::GetZ(); z++) {
                    max_v_z[z] = max_v_z[z] > e.w[z] ? max_v_z[z] : e.w[z];
                }
            }
            for (int z = 0; z < TIC::GetZ(); z++) {
                if (flag == 1)influ[z] += max(max_u_z[z], max_v_z[z]);
                else influ[z] += max_u_z[z] * max_u_z[z];
            }
        }
        double max_influ = 0.0;
        for (int z = 0; z < TIC::GetZ(); z++) {
            if (max_influ < influ[z])
                max_influ = influ[z];
        }
        rt += max_influ;
        fprintf(f, "%d %lg", TIC::GetName(i), rt / TIC::GetCost(TIC::GetName(i)));

        for (int z = 0; z < TIC::GetZ(); z++)
            fprintf(f, " %lg", influ[z]);
        fprintf(f, "\n");
    }
    fclose(f);
}

void keyword_USET(char *filename, double ) {
    double timer;
    struct timeval tstart, tend;
    FILE *f = fopen(filename, "w");
    fprintf(f, "%d\n", TIC::GetN());
    map<int, NB> u_nb;
    double *max_in_z = new double[TIC::GetZ()];
    double *max_out_z = new double[TIC::GetZ()];
    gettimeofday(&tstart, NULL);
    for (int k = 0; k < TIC::GetK(); k++) {
        printf("processing keyword %d...\n", k);
        if (k > 3)break;
        vector<tuple> list = TIM::il_[k].list;
        for (int i = 0; i < list.size(); i++) {
            int v = TIC::GetId(list[i].first);

            double score = list[i].second;
            //get the max in edge of u
            for (int z = 0; z < TIC::GetZ(); z++) max_in_z[z] = 0;
            for (int x = 0; x < TIC::GetInNeighbor(v); x++) {
                Edge e = TIC::GetInEdge(v, x);
                for (int z = 0; z < TIC::GetZ(); z++) {
                    max_in_z[z] = max_in_z[z] > e.w[z] ? max_in_z[z] : e.w[z];
                }
            }
            vector<di> rst = TIC::ReverseDijkstra(v, -log(theta));
            // printf("rst %d\n",rst.size() );
            for (int j = 0; j < rst.size(); j++ ) {
                int u = rst[j].second;
                // printf("u %d name %d\n",u,TIC::GetName(u) );
                // printf("u %d name %d size %d\n",u,TIC::GetName(u),u_nb.size() );
                NB nb;
                map<int, NB>::iterator it;
                it = u_nb.find(TIC::GetName(u));
                if (it == u_nb.end()) {
                    nb.max_ub = 0;
                    for (int z = 0; z < TIC::GetZ(); z++)
                        nb.ub.push_back(0.0);
                } else {
                    nb = u_nb[TIC::GetName(u)];
                }
                if (u == v) {
                    // printf("added u %d k %d %lf %lf\n",u,k,TIM::GetKUWeight(TIC::GetName(u),k),score);
                    nb.max_ub += score;
                    u_nb[TIC::GetName(u)] = nb;
                    continue;
                }

                int flag = 0; //if 1,then v is near u
                for (int z = 0; z < TIC::GetZ(); z++) max_out_z[z] = 0;
                for (int x = 0; x < TIC::GetOutNeighbor(u); x++) {
                    Edge e = TIC::GetOutEdge(u, x);
                    if (v == e.v) {
                        flag = 1;
                    }
                    for (int z = 0; z < TIC::GetZ(); z++) {
                        max_out_z[z] = max_out_z[z] > e.w[z] ? max_out_z[z] : e.w[z];
                    }
                }
                for (int z = 0; z < TIC::GetZ(); z++) {
                    if (flag == 1)nb.ub[z] += max(max_in_z[z], max_out_z[z]) * score;
                    else nb.ub[z] += max_in_z[z] * max_out_z[z] * score;
                }
                u_nb[TIC::GetName(u)] = nb;
            }
        }
        fprintf(f, "%d %d", k, u_nb.size() );
        for (map<int, NB>::iterator it = u_nb.begin(); it != u_nb.end(); it++) {
            NB nb = it->second;
            int u = it->first;
            double max_z = 0.0;
            for (int z = 0; z < TIC::GetZ(); z++)
                if (nb.ub[z] > max_z)max_z = nb.ub[z];
            nb.max_ub += max_z;

            fprintf(f, " %d %.3lf", u, nb.max_ub / TIC::GetCost(u) );
            for (int z = 0; z < TIC::GetZ(); z++)
                fprintf(f, " %.3lf", nb.ub[z] );
        }
        fprintf(f, "\n");
        u_nb.clear();
    }
    gettimeofday(&tend, NULL);
    timer = tend.tv_sec - tstart.tv_sec;
    printf("time %lf\n", timer );
    fclose(f);
}

void USET3(double theta, char *filename, int jobId) {
    double timer;
    struct timeval tstart, tend;
    int span = 1000;
    int length = TIC::GetN() > (jobId + 1) * 1000 ? 1000 : TIC::GetN() - jobId * 1000;
    FILE *f = fopen(filename, "w");
    fprintf(f, "%d\n", length) ;
    map<int, NB> k_nb;
    double *max_u_z = new double[TIC::GetZ()];
    double *max_v_z = new double[TIC::GetZ()];
    double *influ = new double[TIC::GetZ()];
    gettimeofday(&tstart, NULL);
    for (int i = 1000 * jobId; i < 1000 * jobId + length; i++) {
        if (i % 100 == 0) printf("processing vertex %d...\n", i);
        for (int z = 0; z < TIC::GetZ(); z++) max_u_z[z] = 0;
        for (int x = 0; x < TIC::GetOutNeighbor(i); x++) {
            Edge e = TIC::GetOutEdge(i, x);
            for (int z = 0; z < TIC::GetZ(); z++) {
                max_u_z[z] = max_u_z[z] > e.w[z] ? max_u_z[z] : e.w[z];
            }
        }

        vector<di> rst = TIC::Dijkstra(i, -log(theta));
        for (unsigned int j = 0; j < rst.size(); j++) {
            int v = rst[j].second;
            int flag = 0; //if 1,then v is near u
            if (v != i) {
                for (int z = 0; z < TIC::GetZ(); z++) max_v_z[z] = 0;
                for (int z = 0; z < TIC::GetZ(); z++) influ[z] = 0;
                for (int x = 0; x < TIC::GetInNeighbor(v); x++) {
                    Edge e = TIC::GetInEdge(v, x);
                    if (i == e.u)flag = 1;
                    for (int z = 0; z < TIC::GetZ(); z++) {
                        max_v_z[z] = max_v_z[z] > e.w[z] ? max_v_z[z] : e.w[z];
                    }
                }

                for (int z = 0; z < TIC::GetZ(); z++) {
                    if (flag == 1)influ[z] = max(max_u_z[z], max_v_z[z]);
                    else influ[z] = max_u_z[z] * max_v_z[z];
                }
            } else {
                for (int z = 0; z < TIC::GetZ(); z++) influ[z] = 0;
            }
            vector<int>keywords = TIC::GetInterest(TIC::GetName(v));
            for (int k = 0; k < keywords.size(); k++) {
                int keyword = keywords[k];
                double score = TIM::GetKUWeight(TIC::GetName(v), keyword);
                NB nb;
                if (k_nb.count(keyword) == 0) {
                    nb.max_ub = 0.0;
                    for (int z = 0; z < TIC::GetZ(); z++)
                        nb.ub.push_back(0.0);
                } else
                    nb = k_nb[keyword];

                for (int z = 0; z < TIC::GetZ(); z++)
                    nb.ub[z] += influ[z] * score;
                if (v == i)nb.max_ub += score;
                k_nb[keyword] = nb;
            }
        }
        for (map<int, NB>::iterator it = k_nb.begin(); it != k_nb.end(); it++) {
            NB nb = it->second;
            int keyword = it->first;
            int u = TIC::GetName(i);
            double max_z = 0.0;
            for (int z = 0; z < TIC::GetZ(); z++)
                max_z = max_z > nb.ub[z] ? max_z : nb.ub[z];
            nb.max_ub += max_z;

            fprintf(f, "%d %d %.3lf", u, keyword, nb.max_ub / TIC::GetCost(u) );
            for (int z = 0; z < TIC::GetZ(); z++)
                fprintf(f, " %.3lf", nb.ub[z] );
            fprintf(f, "\n");
        }
        k_nb.clear();
    }
    fclose(f);
    gettimeofday(&tend, NULL);
    timer = tend.tv_sec - tstart.tv_sec;
    printf("time %lf\n", timer );
}


int main(int argc, char *argv[]) {
    string cmd;

    cmd = "-nb";
    if (!cmd.compare(argv[1])) {
        TIC::Build2TICFromFile(argv[2]);
        TIC::BuildCost(argv[3]);
        int bound1 = 1000;
        sscanf(argv[4], "%d", &bound1);
        USET(1.0 / bound1, argv[5]);
    }

    cmd = "-knb";
    if (!cmd.compare(argv[1])) {
        TIC::Build2TICFromFile(argv[2]);
        TIC::BuildCost(argv[3]);
        TIC::BuildKeywordTopic(argv[4]);
        TIC::BuildUserInterest(argv[5]);
        TIM::BuildInvertedList(argv[6]);
        int bound1 = 1000;
        sscanf(argv[7], "%d", &bound1);
        // USET2(1.0/bound1,argv[8]);
        USET3(1.0 / bound1, argv[8], 0);
    }

}
