from math import log,exp
import pickle
from random import uniform,randint
import networkx as nx
import time
import numpy as np
# import matplotlib.pyplot as plt  

# path="../data/toy/"
# dataset="toy"
# K=4
# N=8
# Z=2

path="../data/acm/"
dataset="acm"
K=48685
N=1018712
Z=9



def output(filename,data):
    with open(path+"tmp/"+filename,'wb') as fw:
        pickle.dump(data,fw)
def input(filename):
    with open(path+"tmp/"+filename,'r') as fr:
        return pickle.load(fr)
class NB(object):
    """docstring for NB"""
    def __init__(self, max_ub=0.0,vec=[0.0]*Z):
        super(NB, self).__init__()
        self.max_ub = max_ub
        self.vec = vec

def mean_var(data):
    array=np.array(data)
    sum1=array.sum()
    array2=array*array
    sum2=array2.sum()
    N=len(array)
    mean=round(sum1/N,5)
    var=round(sum2/N-mean**2,5)
    print(N,mean,var)             

class Offline(object):
    """docstring for Offline"""
    def __init__(self):
        super(Offline, self).__init__()

    def GenerateCost(self,max_):
        with open("../data/acm/acm.user.interest",'r') as fr:
            N=int(fr.readline().strip())
        with open("../data/acm/acm.cost",'w') as fw:
            fw.write(str(N)+'\n')
            for uid in range(N+1):
                fw.write(str(uid)+' '+str(round(uniform(1,max_),2))+'\n')
    
    def GetInvertedList(self):
        keyword_user={}
        user_keyword={}
        with open(path+dataset+".user.interest",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                if idx %5000==0:print idx
                line=map(int,line.strip().split())
                uid=line[0]
                length=line[1]
                if length<5:continue
                keyword_cnt=user_keyword.get(uid,{})
                for i in range(length):
                    kid=line[2*i+2]
                    cnt=line[2*i+3]
                    if cnt<=1:continue
                    a=keyword_user.get(kid,set())
                    a.add(uid)
                    keyword_user[kid]=a
                    keyword_cnt[kid]=keyword_cnt.get(kid,0)+cnt
                user_keyword[uid]=keyword_cnt
        print "start writing"

        with open(path+dataset+".il5.1",'w') as fw:
            fw.write(str(K)+'\n')
            for k in range(K):
                if k%1000==0:print k
                users=keyword_user.get(k,set())
                fw.write(str(k))
                if len(users)==0:
                    fw.write(" 0 0\n")
                    continue
                idf=log(float(N)/len(users))
                fw.write(" "+str(round(idf,3))+" "+str(len(users)))
                users=list(users)
                users.sort()
                for u in users:
                    keyword_cnt=user_keyword.get(u,{})
                    sum_=0.0
                    for k1 in keyword_cnt:
                        sum_+=keyword_cnt[k1]
                    tf=keyword_cnt[k]/sum_
                    fw.write(" "+str(u)+" "+str(round(tf*idf,3)))
                fw.write('\n')

    def GetKNB(self):
        time1=0
        time2=0
        time3=0
        G=nx.DiGraph()
        time1=time.time()
        with open(path+dataset+".edge",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                line=line.strip().split('\t')
                u=int(line[0])
                v=int(line[1])
                vec=map(float,line[2:])
                G.add_edge(u,v,vec=vec,out_weight=-log(max(vec)))
        print "load in time %lf"%(time.time()-time1)
        time1=time.time()
        for u in G.nodes():
            max_in=[0.0]*Z
            max_out=[0.0]*Z
            for pre in G.predecessors(u):
                max_in=map(lambda (a,b):max(a,b),zip(max_in,G.edges[pre,u]["vec"]))
            for suc in G.successors(u):
                max_out=map(lambda (a,b):max(a,b),zip(max_out,G.edges[u,suc]["vec"]))
            G.nodes[u]["in"]=max_in
            G.nodes[u]["out"]=max_out
        print "cal max in/out time %lf"%(time.time()-time1)

        time1=time.time()
        for u,v in G.edges():
            if not G.has_edge(v,u):continue
            G[v][u]["in_weight"]=G[u][v]["out_weight"]
        print "cal reverse time %lf"%(time.time()-time1)
        il={}
        time1=time.time()
        with open(path+dataset+".il",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                # if idx %1000==0:print idx
                line=line.strip().split()
                kid=int(line[0])
                idf=float(line[1])
                length=int(line[2])
                user_weight=il.get(kid,{})
                for i in range(length):
                    uid=int(line[2*i+3])
                    weight=float(line[2*i+4])
                    user_weight[uid]=weight
                il[kid]=user_weight
        print "user inverted list %lf"%(time.time()-time1)

        user_cost={}
        with open(path+dataset+".cost",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                line=line.strip().split()
                uid=int(line[0])
                cost=float(line[1])
                user_cost[uid]=cost
        time1=0
        start3=time.time()
        with open(path+dataset+".knb",'w') as fw:
            fw.write(str(K)+" "+str(Z)+'\n')
            for kid in range(K):
                if kid %10==0:print kid
                if kid>2:break
                u_nb={}
                user_weight=il[kid]
                print len(user_weight)
                for v in user_weight.keys():
                    max_in=G.nodes[v]["in"]
                    start1=time.time()

                    nodes=nx.single_source_dijkstra(G,source=v,weight="in_weight",cutoff=-log(1e-3))[0].keys()
                    time1+=time.time()-start1
                    start2=time.time()
                    for u in nodes:
                        if u==v:
                            nb = u_nb.get(u,NB())
                            nb.max_ub=nb.max_ub+user_weight[v]
                            u_nb[u]=nb
                            continue
                        max_out=G.nodes[u]["out"]
                        res=[0.0]*9
                        if G.has_edge(u,v):
                            res=map(lambda (a,b):max(a,b),zip(max_in,max_out))
                        else:
                            res=map(lambda (a,b):a*b,zip(max_in,max_out))
                        res=map(lambda a:a*user_weight[v],res)
                        nb = u_nb.get(u,NB())
                        nb.vec=map(lambda (a,b):round(a+b,3),zip(nb.vec,res))
                        nb.max_ub=nb.max_ub+max(res)
                        u_nb[u]=nb
                    time2+=time.time()-start2
                fw.write(str(kid)+" "+str(len(u_nb)))
                for u in u_nb:
                    fw.write(" "+str(u)+" "+str(round(u_nb[u].max_ub/user_cost[u],3))+" "+" ".join(map(str,u_nb[u].vec)))
                fw.write('\n')
            time3=time.time()-start3
            print "dijkstra",time1,"other",time2,"total",time3

    def GetNB(self):
        time1=0
        time2=0
        time3=0
        G=nx.DiGraph()
        time1=time.time()
        with open(path+dataset+".edge",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                line=line.strip().split('\t')
                u=int(line[0])
                v=int(line[1])
                vec=map(float,line[2:])
                G.add_edge(u,v,vec=vec,out_weight=-log(max(vec)))
        print "load in time %lf"%(time.time()-time1)
        time1=time.time()
        for u in G.nodes():
            max_in=[0.0]*Z
            max_out=[0.0]*Z
            for pre in G.predecessors(u):
                max_in=map(lambda (a,b):max(a,b),zip(max_in,G.edges[pre,u]["vec"]))
            for suc in G.successors(u):
                max_out=map(lambda (a,b):max(a,b),zip(max_out,G.edges[u,suc]["vec"]))
            G.nodes[u]["in"]=max_in
            G.nodes[u]["out"]=max_out
        print "cal max in/out time %lf"%(time.time()-time1)
        
        user_cost={}
        with open(path+dataset+".cost",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                line=line.strip().split()
                uid=int(line[0])
                cost=float(line[1])
                user_cost[uid]=cost
        time1=time.time()
        with open(path+dataset+".nb",'w') as fw:
            fw.write(str(N)+" "+str(Z)+'\n')
            for u in range(N):
                if u %10000==0:print u
                nb=NB()
                max_out=G.nodes[u]["out"]
                nodes=nx.single_source_dijkstra(G,source=u,weight="out_weight",cutoff=-log(1e-3))[0].keys()
                for v in nodes:
                    if u==v:
                        nb.max_ub+=1
                        continue
                    max_in=G.nodes[v]["in"]
                    temp=[0.0]*Z
                    if G.has_edge(u,v):
                        temp=map(lambda (a,b):max(a,b),zip(max_in,max_out))
                    else:
                        temp=map(lambda (a,b):a*b,zip(max_in,max_out))
                    nb.vec=map(lambda (a,b):round(a+b,3),zip(nb.vec,temp))
                fw.write(str(u)+" "+str(round(nb.max_ub/user_cost[u],3))+" "+" ".join(map(str,nb.vec))+"\n")
        print "cal  time %lf"%(time.time()-time1)
        
    def GetQuery(self,k=3,foreach=True,length=2,ratio=0.3):
        keyword_topic={}
        topic_keyword={}
        keyword_cnt={}
        with open(path+dataset+".user.interest",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                # if idx %5000==0:print idx
                line=map(int,line.strip().split())
                uid=line[0]
                length2=line[1]
                for i in range(length2):
                    kid=line[2*i+2]
                    cnt=line[2*i+3]
                    keyword_cnt[kid]=keyword_cnt.get(kid,0)+cnt
        keyword_cnt= sorted(keyword_cnt.iteritems(), key=lambda d:d[1], reverse = True)
        keyword_set=set()
        for i in range(len(keyword_cnt)):
            if i>len(keyword_cnt)*ratio:break
            keyword_set.add(keyword_cnt[i][0])
        # print keyword_set
        with open(path+dataset+'.keyword','r') as fr:
            fr.readline()
            fr.readline()
            for line in fr:
                line=line.strip().split()
                kid=int(line[0])
                if kid not in keyword_set:continue
                vec=map(float,line[1:])
                keyword_topic[kid]=vec
                idx=vec.index(max(vec))
                t=topic_keyword.get(idx,[])
                t.append(kid)
                topic_keyword[idx]=t

        with open(path+dataset+'.query','w') as fw:
            cnt=0
            for topic in topic_keyword:
                if len(topic_keyword[topic])<length:continue
                cnt+=1
            fw.write(str(k*cnt)+'\n')
            if not foreach:return
            for topic in topic_keyword:
                t=topic_keyword[topic]
                if len(t)<length:continue
                for i in range(k):
                    cnt=0
                    fw.write(str(length))
                    used=set()
                    while(cnt<length):
                        a=t[randint(0,len(t)-1)]
                        if a in used:continue
                        used.add(a)
                        fw.write(' '+str(a))
                        cnt+=1
                    fw.write('\n')




    # def dijkstra_v_u(self,v,v_u,edge,theta=0.5):
    #     theta=-log(theta)
    #     used_=set()
    #     d={v:0}
    #     while True:
    #         top=self.extract(d)[0]
    #         # if v==2:print top,"-"
    #         if top==-1:break
    #         dist=d[top]
    #         if dist>=theta:break
    #         used_.add(top)
    #         del d[top]
    #         for u in v_u.get(top,set()):
    #             # if v==2:print u,top,max(edge[u][top]),-log(max(edge[u][top]))
    #             if u in used_:continue
    #             w=-log(max(edge[u][top]))
    #             if u not in d or (dist+w)<d[u]:
    #                 d[u]=dist+w
    #     return used_
    # def dijkstra_u_v(self,u,u_v,edge,theta=0.5):
    #     theta=-log(theta)
    #     used_=set()
    #     d={u:0}
    #     while len(d)!=0:
    #         top=min(d,key=d.get)
    #         dist=d[top]
    #         if dist>=theta:break
    #         used_[top]=dist
    #         del d[top]
    #         for v in u_v.get(top,set()):
    #             if v in used_:continue
    #             w=-log(max(edge[top][v]))
    #             if v not in d or (dist+w)<d[v]:
    #                 d[v]=dist+w
    #     return used_
    
    # def extract(self,d):
    #     idx=-1
    #     minimum=float('inf')
    #     for i in d:
    #         if d[i]<minimum:
    #             minimum=d[i]
    #             idx=i
        # return idx,minimum
    # def dijkstra_u_v2(self,u,G,theta=1e-3):
    #     start=time.time()
    #     theta=-log(theta)
    #     used_={}
    #     d={u:0}
    #     while True:
    #         top,dist=self.extract(d)
    #         if top==-1:break
    #         if dist>=theta:break
    #         used_[top]=exp(-dist)
    #         del d[top]
    #         for v in G.successors(top):
    #             if v in used_:continue
    #             w=-log(max(G[top][v]["vec"]))
    #             if v not in d or (dist+w)<d[v]:
    #                 d[v]=dist+w
    #     global dijkstra_time
    #     dijkstra_time+=time.time()-start
    #     return used_
    
    def GetDijkstra(self):
        time1=0
        time2=0
        time3=0
        G=nx.DiGraph()
        #cal max in max out to reduce the heavy load in computation
        time1=time.time()
        with open(path+dataset+".edge",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                # if idx%1000000==0:print idx
                line=line.strip().split('\t')
                u=int(line[0])
                v=int(line[1])
                vec=map(float,line[2:])
                G.add_edge(u,v,vec=vec,out_weight=-log(max(vec)))
        print "load in time %lf"%(time.time()-time1)
        time1=time.time()
        for u in G.nodes():
            max_in=[0.0]*Z
            max_out=[0.0]*Z
            for pre in G.predecessors(u):
                max_in=map(lambda (a,b):max(a,b),zip(max_in,G.edges[pre,u]["vec"]))
            for suc in G.successors(u):
                max_out=map(lambda (a,b):max(a,b),zip(max_out,G.edges[u,suc]["vec"]))
            G.nodes[u]["in"]=max_in
            G.nodes[u]["out"]=max_out
        print "cal max in time %lf"%(time.time()-time1)

        time1=time.time()
        keyword_user={}
        user_keyword={}
        with open(path+dataset+".user.interest",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                # if idx %5000==0:print idx
                line=map(int,line.strip().split())
                uid=line[0]
                length=line[1]
                keyword_cnt=user_keyword.get(uid,{})
                for i in range(length):
                    kid=line[2*i+2]
                    cnt=line[2*i+3]
                    a=keyword_user.get(kid,set())
                    a.add(uid)
                    keyword_user[kid]=a
                    keyword_cnt[kid]=keyword_cnt.get(kid,0)+cnt
                user_keyword[uid]=keyword_cnt
        print "keyword load in in time %lf"%(time.time()-time1)
        
        fw=open(path+dataset+".dijkstra",'w')
        time1=0
        start3=time.time()
        for u in range(N):
            if u>10:break
            if not G.has_node(u):continue
            if u %10000==0:print "%lf"%(float(u)/N*100),"%"
            kid_nb={}
            for k in user_keyword.get(v,{}):
                nb=NB()
                nb.max_ub+=1
                kid_nb[k]=nb
            max_out=G.nodes[u]["out"]
            start1=time.time()
            nodes=nx.single_source_dijkstra(G,u,cutoff=-log(1e-3))[0].keys()
            # nodes=self.dijkstra_u_v2(u,G,theta=1e-3)
            time1+=time.time()-start1
            start2=time.time()
            for v in nodes:
                if u==v:continue
                max_in=G.nodes[v]["in"]
                if G.has_edge(u,v):
                    temp=map(lambda (a,b):max(a,b),zip(max_in,max_out))
                else:
                    temp=map(lambda (a,b):a*b,zip(max_in,max_out))
                for k in user_keyword.get(v,{}):
                    nb=kid_nb.get(k,NB())
                    nb.vec=map(lambda (a,b):(a+b),zip(nb.vec,temp))
                    kid_nb[k]=nb
            time2+=time.time()-start2
            for k in kid_nb:
                nb=kid_nb[k]
                nb.max_ub+=max(nb.vec)
                nb.vec=map(lambda a:round(a,3),nb.vec)
                with open(path+dataset+"/tmp/"+str(k)+".dijkstra",'a') as fw: 
                    fw.write(str(u)+" "+str(nb.max_ub)+" "+" ".join(map(str,nb.vec))+"\n")
        time3=time.time()-start3
        print "dijkstra",time1,"other",time2,"total",time3

   
    def PlotTF(self):
        tf=[]
        tfidf=[]
        with open(path+dataset+".il5.1",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                if idx %5000==0:print idx
                line=line.strip().split()
                idf=float(line[1])
                length=int(line[2])
                for i in range(length):
                    weight=float(line[2*i+4])
                    # if weight/idf>0.5:print weight/idf,length
                    tf.append(weight/idf)
                    tfidf.append(weight)
                if idx>0000:break
        mean_var(tf)
        mean_var(tfidf)
        # sum2
        # print(round(float(sum(tf)) / max(len(tf), 1), 3))
        # print
        # weights=np.ones_like(tf)/float(len(tf))
        # # tf=np.array(tf)
        # # print tf
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # # print tf
        # ax.hist(tf, 100, facecolor='blue',weights=weights,cumulative=True) 
        # # ax.hist(tf, 100, facecolor='blue',weights=weights) 
        # plt.show()  

    def KeywordStats(self):
        u_cnt=0.0
        keyword_cnt=0.0
        unique_keyword_cnt=0.0
        with open(path+dataset+".user.interest",'r') as fr:
            fr.readline()
            for idx,line in enumerate(fr):
                if idx %100000==0:print idx
                line=map(int,line.strip().split())
                uid=line[0]
                length=line[1]
                
                if length<5:continue
                u_cnt+=1
                unique_keyword_cnt+=length
                for i in range(length):
                    kid=line[2*i+2]
                    cnt=line[2*i+3]
                    if cnt<=3:unique_keyword_cnt-=1
                    else:keyword_cnt+=cnt
        print(u_cnt,keyword_cnt,unique_keyword_cnt)
        print(keyword_cnt/u_cnt,unique_keyword_cnt/u_cnt)


    def process_user_interest(self,length_min=5,keyword_cnt_min=1):
        fw=open(path+dataset+'user.interest.2','w')
        with open(path+dataset+".user.interest",'r') as fr:
            fw.write(fr.readline())
            for idx,line in enumerate(fr):
                if idx %5000==0:print idx
                line=map(int,line.strip().split())
                uid=line[0]
                length=line[1]
                keywords=[]
                for i in range(length):
                    kid=line[2*i+2]
                    cnt=line[2*i+3]
                    if cnt<=keyword_cnt_min:continue
                    keywords.append((kid,cnt))
                if len(keywords)<=length_min:continue
                fw.write(str(uid)+" "+str(len(keywords)))
                for kid,cnt in keywords:
                    fw.write(" "+str(kid)+" "+str(cnt))
                fw.write("\n")
        fw.close()


if __name__ == '__main__':
    offline=Offline()
    offline.PlotTF()
    # offline.GetInvertedList()
    # offline.KeywordStats()
    # offline.GetKNB()
    # offline.GetNB()
    # offline.GenerateCost(3)
    # offline.GetQuery()
    # offline.GetDijkstra()






        