#include "stdio.h"
#include "sys/time.h"
#include <string>
#include <math.h>

#include "tim.h"

void PrintS(S s) {
    printf("S = \n");
    for (unsigned int i = 0; i < s.ids.size(); i++) {
        printf("%d\t%lg\n", s.ids[i], s.minflu[i]);
    }
}

int main(int argc, char *argv[]) {
    double timer;

    struct timeval tstart,tend;

    string cmd;
    // cmd = "-g";//Greedy
    // if(!cmd.compare(argv[1])){
    //     TIM::LoadTIC(argv[2],argv[3]);
    //     int bound1, bound2;
    //     TIM::Offline(argv[4]);
    //     sscanf(argv[5], "%d", &bound1);
    //     sscanf(argv[6], "%d", &bound2);
    //     TIM::theta_ = 1.0/bound1;
    //     TIM::theta2_ = 1.0/bound2;
    //     FILE *f = fopen(argv[7], "r");
    //     int nq = 0, Z = 0;
    //     double B;
    //     sscanf(argv[8], "%lf", &B);

    //     vector<double> item;
    //     fscanf(f, "%d %d", &nq, &Z);
    //     item.resize(Z);
    //     for (int i = 0; i < nq; i++) {
    //         gettimeofday(&tstart,NULL);
    //         for (int z = 0; z < Z; z++) {
    //             fscanf(f, "%lg", &(item[z]));
    //         }

    //         Q q;
    //         q.item = item;
    //         q.B = B;
    //         S s;

    //         s = TIM::Greedy(q);

    //         gettimeofday(&tend,NULL);
    //         timer=1000000*(tend.tv_sec-tstart.tv_sec)+tend.tv_usec-tstart.tv_usec;

    //         vector<int> temp;
    //         for(int i=0;i<s.ids.size();i++)
    //             temp.push_back(TIC::GetId(s.ids[i]));
    //         double spread = TIC::Spread(temp, q.item);
            
    //         double rt = 0;
    //         for (unsigned int i = 0; i < s.ids.size(); i++) {
    //             rt += s.minflu[i];
    //         }

    //         #ifdef VERBOSE
    //         printf("spread : %lg, time: %lg S size:%d\n", spread, timer/CLOCKS_PER_SEC,s.ids.size());

    //         #endif
    //         #ifdef EXP
    //         printf("%lg %lg %lg\n", timer/CLOCKS_PER_SEC, spread, B);
    //         #endif
            
    //     }
    //     fclose(f);
    // }
    cmd = "-b";   //Besteffort
    if (!cmd.compare(argv[1])){       
        TIC::Build2TICFromFile(argv[2]);
        TIC::BuildCost(argv[3]);
        // printf("here\n");
        TIC::BuildUserInterest(argv[4]);
        // printf("here\n");
        TIC::BuildKeywordTopic(argv[5]);
        // printf("here\n");
        int bound1, bound2;
        sscanf(argv[10], "%d", &bound1);
        sscanf(argv[11], "%d", &bound2);
        TIM::theta_ = 1.0/bound1;
        TIM::theta2_ = 1.0/bound2;
        // printf("here\n");
        TIM::BuildInvertedList(argv[6]);
        // printf("here\n");
        TIM::BuildKeywordNeighborhoodBound(argv[7]);
        printf("end loading\n"); 
        TIM::Init();       
        FILE *f = fopen(argv[8], "r");
        int nq = 0, Z = 0;
        double B;
        sscanf(argv[9], "%lf", &B);      
        fscanf(f, "%d", &nq);
        for (int i = 0; i < nq; i++) {
            vector<int> keywords;
            gettimeofday(&tstart,NULL);
            int num=0;
            fscanf(f, "%d", &num);
            for(int j=0;j<num;j++){
                int kid;
                fscanf(f,"%d",&kid);
                keywords.push_back(kid);
            }
            Q q;
            q.item=TIM::GetItem(keywords);
            q.B = B;
            q.keywords=keywords;
            S s;
            s = TIM::BestEffort(q);
            printf("result\n");
            gettimeofday(&tend,NULL);
            timer=1000000*(tend.tv_sec-tstart.tv_sec)+tend.tv_usec-tstart.tv_usec;

            vector<int> ids,names;
            for(int i=0;i<s.ids.size();i++){
                ids.push_back(TIC::GetId(s.ids[i]));
                names.push_back(s.ids[i]);
            }
            double spread = TIC::Spread(ids, q.item);
            double score = TIM::WeightedScore(names,q);
            double rt = 0;
            for (unsigned int i = 0; i < s.ids.size(); i++) {
                rt += s.minflu[i];
            }

            // #ifdef VERBOSE
            // printf("spread : %lg, time: %lg S size:%d\n", spread, timer/CLOCKS_PER_SEC,s.ids.size());
            // // PrintS(s);
            // #endif
            #ifdef EXP
            printf("%lg %lg %lg %lg %lg\n", timer/CLOCKS_PER_SEC, rt, B,spread,score);
            #endif
            
        }
        fclose(f);
    }

    cmd = "-n";   //neighborhood
    if (!cmd.compare(argv[1])){       
        TIC::Build2TICFromFile(argv[2]);
        TIC::BuildCost(argv[3]);
        TIC::BuildUserInterest(argv[4]);
        TIC::BuildKeywordTopic(argv[5]);
        int bound1, bound2;
        sscanf(argv[10], "%d", &bound1);
        sscanf(argv[11], "%d", &bound2);
        TIM::theta_ = 1.0/bound1;
        TIM::theta2_ = 1.0/bound2;
        TIM::BuildInvertedList(argv[6]);
        TIM::BuildNeighborhoodBound(argv[7]);
        printf("end loading\n"); 
        TIM::Init();       
        FILE *f = fopen(argv[8], "r");
        int nq = 0, Z = 0;
        double B;
        sscanf(argv[9], "%lf", &B);
        fscanf(f, "%d", &nq);
        for (int i = 0; i < nq; i++) {
            vector<int> keywords;
            gettimeofday(&tstart,NULL);
            int num=0;
            fscanf(f, "%d", &num);
            for(int j=0;j<num;j++){
                int kid;
                fscanf(f,"%d",&kid);
                keywords.push_back(kid);
            }
            Q q;
            q.item=TIM::GetItem(keywords);
            q.B = B;
            q.keywords=keywords;
            S s;
            s = TIM::Neighborhood(q);
            printf("result\n");
            gettimeofday(&tend,NULL);
            timer=1000000*(tend.tv_sec-tstart.tv_sec)+tend.tv_usec-tstart.tv_usec;

            vector<int> temp;
            for(int i=0;i<s.ids.size();i++)
                temp.push_back(TIC::GetId(s.ids[i]));
            double spread = TIC::Spread(temp, q.item);
            double score = TIM::WeightedScore(temp,q);
            double rt = 0;
            for (unsigned int i = 0; i < s.ids.size(); i++) {
                rt += s.minflu[i];
            }

            // #ifdef VERBOSE
            // printf("spread : %lg, time: %lg S size:%d\n", spread, timer/CLOCKS_PER_SEC,s.ids.size());
            // // PrintS(s);
            // #endif
            #ifdef EXP
            printf("%lg %lg %lg %lg %lg\n", timer/CLOCKS_PER_SEC, rt, B,spread,score);
            #endif
            
        }
        fclose(f);
    }
    cmd = "-t";
    if(!cmd.compare(argv[1])){
        TIC::BuildKeywordTopic(argv[2]);
        FILE *f = fopen(argv[3], "r");
        int nq = 0;
        vector<int> keywords;
        fscanf(f, "%d", &nq);
        printf("nq %d\n",nq );
        for (int i = 0; i < nq; i++) {
            int num=0;
            fscanf(f, "%d", &num);
            for(int j=0;j<num;j++){
                int kid;
                fscanf(f,"%d",&kid);
                keywords.push_back(kid);
            }
            Q q;
            q.item=TIM::GetItem(keywords);
            keywords.clear();
        }
    }
    // cmd = "-a";  // Approximation.
    // if (!cmd.compare(argv[1])) {
    //     TIM::LoadTIC(argv[2],argv[3]);
    //     int bound1, bound2;
    //     sscanf(argv[7], "%d", &bound1);
    //     sscanf(argv[8], "%d", &bound2);
    //     TIM::theta_ = 1.0/bound1;
    //     TIM::theta2_ = 1.0/bound2;
    //     TIM::InitApprox(argv[4],argv[5],argv[6]);
        
    //     FILE *f = fopen(argv[9], "r");
    //     int nq = 0, Z = 0;
    //     double B,epsilon;
    //     sscanf(argv[10], "%lf", &epsilon);
    //     sscanf(argv[11], "%lf", &B);
    //     // printf("B:%lf\n", B);
    //     vector<double> item;
    //     fscanf(f, "%d %d", &nq, &Z);
    //     item.resize(Z);
    //     for (int i = 0; i < nq; i++) {
    //         //start = clock();
    //         gettimeofday(&tstart,NULL);
    //         for (int z = 0; z < Z; z++) {
    //             fscanf(f, "%lg", &(item[z]));
    //         }
            
    //         Q q = {.item=item, .B=B};
    //         S s = TIM::Approx(q, epsilon);
           


    //         //end = clock();
    //         gettimeofday(&tend,NULL);
    //         //timer = double(end - start);
    //         timer=1000000*(tend.tv_sec-tstart.tv_sec)+tend.tv_usec-tstart.tv_usec;
    //         double spread = TIC::Spread(s.ids, q.item);

    // #ifdef VERBOSE
    //         printf("spread : %lg, time: %lg\n", spread, timer/CLOCKS_PER_SEC);
    // #endif
    // #ifdef EXP
    //     printf("%lg %lg %lg %lg\n", timer/1000000, spread, B, epsilon);
    // #endif
    //     }
    //     fclose(f);
    // }


}
