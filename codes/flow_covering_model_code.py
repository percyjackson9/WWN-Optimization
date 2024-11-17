### this file contains only codes for flow covering model --- must contain all relevant python functions
"""
This code is specific for LA data and a test data-set
"""
import cplex
from cplex.exceptions import CplexError
from cplex.callbacks import LazyConstraintCallback as LCC
from cplex.callbacks import BranchCallback as BC
from cplex.callbacks import NodeCallback as NC
from cplex.callbacks import UserCutCallback as UCC
from cplex.callbacks import MIPInfoCallback as INFO
from cplex.callbacks import Context
import numpy as np
import copy
import collections
import docplex.mp.environment as envi
import docplex.mp.model as dp
from docplex.mp.callbacks.cb_mixin import *

import time
import networkx as nx
from shapely.geometry import Point,LineString
import geopandas
import matplotlib.pyplot as plt
import pickle
import pandas as pd

from pathlib import PurePath
import os
import pickle
import sys

'''
Code description:
This code contains all the functions needed to solve the flow covering model. The code is divided into several parts:
instance name can be either "LA" or "regina".
For "LA" the number of clusters is fixed to 7, and for "regina" the number of clusters is fixed to 10.
'''

curd = PurePath(os.path.dirname(os.path.abspath(__file__)))
pdir = curd.parent
picdir = os.path.join(pdir,"pickle_files")
la_dir = os.path.join(picdir,"la_network")
regina_dir = os.path.join(picdir,"regina_network")
results_dir = os.path.join(pdir,"results")
lp_file_dir = os.path.join(pdir,"lp_files")
lpdir = lp_file_dir
nodefile_dir = lp_file_dir


class create_task:
    LA_CLUSTERS = 7
    REGINA_CLUSTERS = 10
    ALGO_MAPPER = {"M1":"basicnored",
                    "M2":"basic",
                    "M3":"oneshot",
                    "M4":"contextwithlb"}

    LA_ZONE_MAPPER = {0:0,
                      1:2,
                      2:1}      
    
    REGINA_ZONE_MAPPER = {1:1,  
                          2:4,
                          3:5,
                          4:8,
                          5:9}  

    def __init__(self,inst_name: str,task: list):
        self.inst_name=inst_name
        self.task=task
        self._generate_task()
    
    def _generate_task(self):
        self.new_task = [create_task.ALGO_MAPPER[self.task[0]],self.task[1],self.task[2]]
        if self.inst_name=="LA" and len(self.task)==4:
            self.new_task.append(create_task.LA_CLUSTERS)
            new_zones = [create_task.LA_ZONE_MAPPER[i] for i in self.task[3]]
            self.new_task.append(new_zones)
        if self.inst_name=="regina" and len(self.task)==4:
            self.new_task.append(create_task.REGINA_CLUSTERS)
            new_zones = [create_task.REGINA_ZONE_MAPPER[i] for i in self.task[3]]
            self.new_task.append(new_zones)
            
    def get_task(self):
        return self.new_task


class global_params:
    """
    this class stores the solutions
    """
    def __init__(self,id,tl=7200):
        self.id=id
        self.time_limit=tl
        self.sol_ctr=0
        self.solution_queue={}
        self.current_best=0
        self.dual_times=[0,0] ## [count,time]
        self.primal_times=[0,0] ## [count,time]
        self.bender=[0,0] ## [count,time]

    def __str__(self):
        return f"Global parameters for instance specification {self.id}"

    ## this function is used to save solution in results folder
    def save_csv_solution(self,prob,data,task,selected,flt,houses,t,fl):
        self.prob=prob
        infeas = [103, 106, 108, 112]
        num_mnodes=len([i for i in data.G.nodes if data.G.nodes[i]["type"]=="M"])
        num_n=len(data.G.nodes)
        num_e=len(data.G.edges)
        num_h=houses
        if prob.solution.get_status() not in infeas:
            data_to_write = {"ON": [len(data.trunc_nodes)],
                             "OE": [len(data.trunc_edges)],
                             "OM": [len(data.trunc_mnodes)],
                             "OH": [num_h],
                             "N": [num_n],
                             "E": [num_e],
                             "M": [num_mnodes],
                             "H": [num_h],
                             "Time": [t],
                             "Obj": [prob.solution.get_objective_value()],
                             "Gap": [prob.solution.MIP.get_mip_relative_gap()],
                             "Status": [prob.solution.get_status()],
                             "Locs": [selected],
                             "Flows": [flt]}


        elif task[0]=="lazyrej" or task[0]=="contextrej":
            best_bound=prob.solution.MIP.get_best_objective()
            gap=(self.current_best-best_bound)/(self.current_best+1e-20)
            data_to_write = {"ON": [len(data.trunc_nodes)],
                             "OE": [len(data.trunc_edges)],
                             "OM": [len(data.trunc_mnodes)],
                             "OH": [num_h],
                             "N": [num_n],
                             "E": [num_e],
                             "M": [num_mnodes],
                             "H": [num_h],
                             "Time": [t],
                             "Obj": [self.current_best],
                             "Gap": [gap],
                             "Status": [prob.solution.get_status()],
                             "Locs": [gpars.solution_queue[self.sol_ctr][0]],
                             "Flows": [gpars.solution_queue[self.sol_ctr][1]]}

        else:
            data_to_write ={"ON": [len(data.trunc_nodes)],
                             "OE": [len(data.trunc_edges)],
                             "OM": [len(data.trunc_mnodes)],
                             "OH": [num_h],
                             "N": [num_n],
                             "E": [num_e],
                             "M": [num_mnodes],
                             "H": [num_h],
                             "Time": ["INF"],
                             "Obj": ["INF"],
                             "Gap": ["INF"],
                             "Status": [prob.solution.get_status()],
                             "Locs": ["INF"],
                             "Flows": ["INF"]}

        if "lb" in task[0]:
            data_to_write["WST"]=[self.wst]
            data_to_write["WLB"]=[self.newlb]
            data_to_write["WUB"]=[self.newub]


        data_to_write = pd.DataFrame(data_to_write)
        data_to_write.to_csv(os.path.join(results_dir,f"{self.id}.csv"))
        if fl!=[] and task[0]!="lazyrej":
            with open(os.path.join(results_dir,f"{self.id}_fvals.pykl"),"wb") as f:
                pickle.dump(fl,f)

    ## this function notes down run times for warm start
    def save_ws_time(self,t,newlb,newub):
        self.wst=t
        self.newub=newub
        self.newlb=newlb

class ValidInstanceName:
    def __set_name__(self, owner, name):
        self.name=name
    
    def __get__(self,instance,owner):
        return instance.__dict__[self.name]
    
    def __set__(self,instance,value):
        if value not in ["LA","regina"]:
            raise ValueError("Instance name must be either LA or regina")
        instance.__dict__[self.name]=value

class ReadData:
    name = ValidInstanceName()
    def __init__(self,name,base,noc=0,target=[]):
        """
        only accepts name as "LA" or "regina".
        
        """
        self.name=name
        if noc>0:
            if name == "LA":
                assert noc==7, "Number of clusters must be 7 for LA data."
            if name == "regina":
                assert noc==10, "Number of clusters must be 10 for regina data."

            self.clusters=noc
            self.target=target[:]
            if target == []:
                raise Exception("Number of clusters is provided, but target is empty. Please re-run with target as a list of areas.")
            
        self.base=base
        
        if name=="LA":
            print("Getting LA")
            self.get_LA_network()
        if name=="regina":
            print("Getting regina")
            self.get_regina_network()
        
    def generate_solution_shapefile(self,task):
        if len(task)==3:
            fprefix = f"{self.name}_{task[0]}_{task[1]}_{task[2]}"
        else:
            fprefix = f"{self.name}_{task[0]}_{task[1]}_{task[2]}_{task[4]}"
        filename=f"{fprefix}.csv"
        try:
            data=pd.read_csv(os.path.join(results_dir,filename))
        except FileNotFoundError as e:
            print("Must have completed the task first: ",e)
            exit()

                
        locs=data.Locs.values[0].replace("'"," ")
        locs=locs.replace("["," ")
        locs=locs.replace("]"," ")
        locs=locs.replace(","," ")
        LOCS=[int(i) for i in locs.split(" ") if len(i)>0]

        if inst_name=="LA":
            f1 = open(os.path.join(la_dir,"LA_new_cord.pickle"), 'rb')
            node_xy = pickle.load(f1)
            f1.close()
        if inst_name=="regina":
            f1 = open(os.path.join(regina_dir,"regina_all_cord.pickle"), 'rb')
            node_xy = pickle.load(f1)
            f1.close()

        PL=[Point(node_xy[i]) for i in LOCS]
        newdata = geopandas.GeoDataFrame(geometry=PL)
        newdata.to_file(os.path.join(results_dir,f"{fprefix}_solution.shp"))

        if len(task)==3:
            PL=[Point(node_xy[i]) for i in node_xy.keys()]
            newdata = geopandas.GeoDataFrame(geometry=PL)
            newdata.to_file(os.path.join(results_dir,f"{fprefix}_full_network.shp"))
        else:
            PL=[Point(node_xy[i]) for i in node_xy.keys() if self.cluster_select[i]==True]
            newdata = geopandas.GeoDataFrame(geometry=PL)
            newdata.to_file(os.path.join(results_dir,f"{fprefix}_full_network.shp"))
    
        print("Shapefiles saved in results folder")

    ## this function plots a set of nodes.
    def plot_nodexy(self,nodes):
        fig=plt.figure(1)
        for i in nodes:
            plt.plot(self.node_xy[i][0],self.node_xy[i][1],markersize=1,color='r',marker="o")
        plt.title(f"nodes for {self.name}")
        plt.axis("scaled")
        plt.show()

    ## this function stores selected nodes for sensors.
    def store_solution(self,sel):
        self.selected=sel[:]

    ## this function plots the network with selected nodes.    
    def plot_network(self,nodes,mnodes, edges, node_xy,selected=None, target= [], houses=[],text=None, to_remove=None):
        if selected is None:
            selected = []
        for i in edges:
            plt.plot([node_xy[i[0]][0], node_xy[i[1]][0]], [node_xy[i[0]][1], node_xy[i[1]][1]], "k")
        if to_remove != [] and to_remove != None:
            for i in to_remove:
                plt.plot([node_xy[i[0]][0], node_xy[i[1]][0]], [node_xy[i[0]][1], node_xy[i[1]][1]], "r")

        for i in nodes:
            if i in target:
                plt.plot(node_xy[i][0], node_xy[i][1], "ro")
                plt.text(node_xy[i][0], node_xy[i][1], s=f"{i}")
            elif i in houses:
                plt.plot(node_xy[i][0], node_xy[i][1], "go")
                plt.text(node_xy[i][0], node_xy[i][1], s=f"{i}")
            elif i in mnodes:
                plt.plot(node_xy[i][0], node_xy[i][1], "bo")
                plt.text(node_xy[i][0], node_xy[i][1], s=f"{i}")
            else:
                plt.text(node_xy[i][0], node_xy[i][1], s=f"{i}")

        if selected != None and selected != []:
            for i in selected:
                plt.plot(node_xy[i][0], node_xy[i][1], "co", markersize=10)
            plt.title(f"Network with solution, {text}")
        else:
            plt.title(f"Original network : {text}")
        plt.axis('scaled')
        plt.show()

    def plot_network_based_on_boundary(self,LB,selected=None, target= [], houses=[],text=None, to_remove=None):
        if selected is None:
            selected = []
        boundb=[]
        bounda=[]
        for i in self.trunc_edges:
            plt.plot([self.node_xy[i[0]][0], self.node_xy[i[1]][0]], [self.node_xy[i[0]][1], self.node_xy[i[1]][1]], "k")
            if self.tr_mpots[i[0]]>=LB and self.tr_mpots[i[1]]<=LB:
                boundb.append(i[0])
                bounda.append(i[1])

        if to_remove != [] and to_remove != None:
            for i in to_remove:
                plt.plot([self.node_xy[i[0]][0], self.node_xy[i[1]][0]], [self.node_xy[i[0]][1], self.node_xy[i[1]][1]], "r")

        if selected != None and selected != []:
            for i in selected:
                plt.plot(self.node_xy[i][0], self.node_xy[i][1], "co", markersize=10)
            plt.title(f"Network with solution, {text}")
        else:
            plt.title(f"Original network : {text}")

        for i in self.trunc_nodes:
            if i in boundb:
                plt.plot(self.node_xy[i][0], self.node_xy[i][1], "go")
                plt.text(self.node_xy[i][0], self.node_xy[i][1], s=f"({i},{self.tr_mpots[i]})")
            elif i in bounda:
                plt.plot(self.node_xy[i][0], self.node_xy[i][1], "ro")
                plt.text(self.node_xy[i][0], self.node_xy[i][1], s=f"({i},{self.tr_mpots[i]})")
            elif i in houses:
                plt.plot(self.node_xy[i][0], self.node_xy[i][1], "go")
                plt.text(self.node_xy[i][0], self.node_xy[i][1], s=f"{i}")
            elif i in self.trunc_mnodes:
                plt.plot(self.node_xy[i][0], self.node_xy[i][1], "bo")
                plt.text(self.node_xy[i][0], self.node_xy[i][1], s=f"({i},{self.tr_mpots[i]})")
            else:
                plt.text(self.node_xy[i][0], self.node_xy[i][1], s=f"{i}")


        plt.axis('scaled')
        plt.show()

    ## exclusively used for getting network features.
    def dynamic_way(self,nodes, adj, dist, edges):
        '''
        adj: adjacency list of a directed network, where direction is opposite to the flow of sewage.
        pred: predecessor of a node, if a node has no predecessor, then it is a source node.
        unc_set: set of nodes that get disconnected from the base node, if egde (i,j) is removed.
        '''
        unc_set = collections.defaultdict(list)
        pred = collections.defaultdict(int)
        print("Started preprocessing")
        for i in nodes:
            if adj[i] != []:
                for j in adj[i]:
                    pred[j] = i
        for i in nodes:
            if i not in pred.keys():
                pred[i] = -1
        for i in nodes:
            if adj[i] == []:
                unc_set[pred[i], i].append(i)
        print("ended preprocessing")
        check = collections.defaultdict(int)
        maxd = max(dist.values())
        for i in list(range(maxd - 1, 0, -1)):
            a = [j for j in edges if dist[j[1]] == i]
            for j in a:
                # t=time.time()
                for k in adj[j[1]]:
                    unc_set[j] += unc_set[j[1], k]
                    check[j[1], k] += 1

        return unc_set, check

    @staticmethod
    def get_directed_rooted_tree(nodes,adj,root):
        pred=collections.defaultdict(int)
        visited=collections.defaultdict(int)
        nn=len(nodes)
        pred[root]=-1
        s=[root]
        dist=collections.defaultdict(lambda:nn)
        dist[root]=0
        while len(s)>0:
            e=s[0]
            visited[e] = 1
            for i in adj[e]:
                if visited[i]==0:
                    dist[i]=dist[e]+1
                    s.append(i)
                    pred[i]=e
                    visited[i]=1
            s.remove(e)

        for i in nodes:
            if adj[i]!=[]:
                adj[i]=[k for k in adj[i] if pred[k]==i]
        tbr=[]
        for i in nodes:
            if dist[i]>=nn:
                adj.pop(i)
                tbr.append(i)
        nodes=list(i for i in nodes if dist[i]<nn)
        for i in tbr:
            dist.pop(i)
        edges=[]
        for i in nodes:
            for j in adj[i]:
                edges.append((i,j))
        
        return nodes,adj,dist,edges

    ## this function is used to truncate whole networks.
    @staticmethod   
    def truncate_network_using_edges_fast(nodes, adj, mnodes, base, edges, npots, weights=None):
        t=time.time()
        """
        Supposed to be a fast implementation
        :param nodes:
        :param adj:
        :param mnodes:
        :param base:
        :param edges:
        :param npots:
        :param weights:
        :return:
        """
        ## npots here includes manhole potentials (node potentials)
        if weights == None:
            weights = {}
            for i in nodes:
                if adj[i] == []:
                    weights[i] = 1
                else:
                    weights[i] = 0

        to_remove = {i: True for i in edges}

        is_mnode = collections.defaultdict(bool)
        for i in mnodes:
            is_mnode[i] = True

        for i in edges:
            if is_mnode[i[0]] and is_mnode[i[1]]:
                if npots[i[0]] != npots[i[1]] and npots[i[1]] > 0:
                    to_remove[i] = False
            if is_mnode[i[1]] and not is_mnode[i[0]]:
                to_remove[i] = False
            if i[0] == base:
                to_remove[i] = False


        if len(to_remove) > 0:
            edges = [i for i in edges if to_remove[i]]

        # edges=list(set(edges).difference(set(to_remove)))
        #print("Time to get list of edges: ", round(time.time() - t, 20))
        print("Starting iteration", len(edges))
        while len(edges) > 0:
            to_change = edges[0]
            weights[to_change[0]] += weights[to_change[1]]
            adj[to_change[0]].remove(to_change[1])
            # print("To change ",to_change)
            # print("Adj ",adj[to_change[1]])
            for i in adj[to_change[1]]:
                to_remove[(to_change[0], i)] = False
                if to_remove[(to_change[1], i)]:
                    edges.remove((to_change[1], i))
                    #to_remove[(to_change[1], i)]=False
                    if is_mnode[to_change[0]] and is_mnode[i]:
                        if npots[to_change[0]] == npots[i] or npots[i] == 0:
                            edges.append((to_change[0], i))
                            to_remove[(to_change[0], i)] = True

                    else:
                        edges.append((to_change[0], i))
                        to_remove[(to_change[0], i)] = True

                else:
                    if is_mnode[to_change[0]] and is_mnode[i]:
                        if npots[to_change[0]] == npots[i] or npots[i]==0:
                            edges.append((to_change[0], i))
                            to_remove[(to_change[0], i)] = True

                adj[to_change[0]].append(i)

            nodes.remove(to_change[1])
            if is_mnode[to_change[1]]:
                mnodes.remove(to_change[1])
            adj.pop(to_change[1])
            weights.pop(to_change[1])
            edges.remove(to_change)

        new_edges = [(i, j) for i in nodes for j in adj[i]]
        print(f"Truncation time {time.time()-t}")

        return weights, nodes, adj, mnodes, new_edges

    def get_los_angeles_clusters(self,noc=7):
        seed=1
        picdir = la_dir
        name = "los_angeles"
        ##### --------------- Regina Data -------------------- #####
        
        f1=open(os.path.join(picdir,"LA_new_cord.pickle"),'rb')
        node_xy=pickle.load(f1)
        f1.close()

        nodes=list(node_xy.keys())

        fname=name+"_"+str(noc)+".pickle"
        points=np.array([node_xy[i] for i in nodes])
        id={i:nodes[i] for i in range(len(nodes))}

        try:
            print("Trying")
            f1 = open(os.path.join(picdir,fname), 'rb')
            cluster_label = pickle.load(f1)
            labels=np.array([cluster_label[i] for i in nodes])
            f1.close()
            return cluster_label,node_xy

        except Exception as e:
            print("Following exception occured: ",e)
            
    def get_regina_clusters(self,noc=10):
        seed=1
        picdir = regina_dir
        name = "regina"
        ##### --------------- Regina Data -------------------- #####
        #"""
    
        f1=open(os.path.join(picdir,"regina_all_cord.pickle"),'rb')
        node_xy=pickle.load(f1)
        f1.close()

        nodes=list(node_xy.keys())

        fname=name+"_"+str(noc)+".pickle"
        points=np.array([node_xy[i] for i in nodes])
        id={i:nodes[i] for i in range(len(nodes))}

        try:
            print("Trying")
            f1 = open(os.path.join(picdir,fname), 'rb')
            cluster_label = pickle.load(f1)
            labels=np.array([cluster_label[i] for i in nodes])
            f1.close()
            return cluster_label,node_xy
        except Exception as e:
            print("Following exception occured: ",e)

    def get_regina_network(self):
        print("Trying to get regina")
        t = time.time()
        f1 = open(os.path.join(regina_dir,"regina_all_cord.pickle"), 'rb')
        f2 = open(os.path.join(regina_dir,"regina_all_adj.pickle"), "rb")
        f3 = open(os.path.join(regina_dir,"regina_all_manholes.pickle"), "rb")
        self.node_xy = pickle.load(f1)
        self.adj = pickle.load(f2)
        self.mnodes = pickle.load(f3)
        f1.close()
        f2.close()
        f3.close()

        self.nodes = list(self.node_xy.keys())

        dedges=[(i,j) for i in self.adj.keys() for j in self.adj[i]]
        edges=set()
        for (i,j) in dedges:
            edges.add((i,j))
            edges.add((j,i))
        edges=list(edges)
        G=nx.Graph()
        G.add_nodes_from(self.nodes)
        G.add_edges_from(dedges)
        self.adj=nx.to_dict_of_lists(G)
        adj=copy.deepcopy(self.adj)
        t=time.time()
        self.nodes, self.adj, dist, self.edges = ReadData.get_directed_rooted_tree(self.nodes, adj, self.base)
        sub_mnodes = list(set(self.nodes).intersection(set(self.mnodes)))
        print(f"Total nodes {len(self.nodes)} total mnodes {len(sub_mnodes)}")
        print(f"DRT in {time.time()-t}")
        t=time.time()
        cut_set, check = self.dynamic_way(self.nodes, self.adj, dist, self.edges)
        print(f"MPOTS in {time.time()-t}")

        self.weights = {}
        for i in self.nodes:
            if self.adj[i] == []:
                self.weights[i] = 1
            else:
                self.weights[i] = 0
        Houses = sum([self.weights[i] for i in self.weights.keys()])

        self.wye_set = {}
        for i in self.edges:
            self.wye_set[i] = [j for j in cut_set[i] if self.weights[j] > 0]

        t = time.time()
        self.mpots = collections.defaultdict(int)
        for i in sub_mnodes:
            if self.adj[i] != []:
                for j in self.adj[i]:
                    self.mpots[i] += len(self.wye_set[(i, j)])
            else:
                self.mpots[i] = 0

        self.npots = collections.defaultdict(int)
        for i in self.nodes:
            if self.adj[i] != []:
                for j in self.adj[i]:
                    self.npots[i] += len(self.wye_set[i, j])
            else:
                self.npots[i] = 0

        self.mnodes=sub_mnodes[:]

    def get_LA_network(self):
        print("Reading LA network")
        f1 = open(os.path.join(la_dir,"LA_new_cord.pickle"), 'rb')
        f2 = open(os.path.join(la_dir,"LA_tree_adj.pickle"), "rb")
        f3 = open(os.path.join(la_dir,"LA_new_mnodes.pickle"), "rb")
        f4 = open(os.path.join(la_dir,"LA_new_wyes.pickle"), "rb")
        f5 = open(os.path.join(la_dir,"LA_tree_pred.pickle"), "rb")
        f6 = open(os.path.join(la_dir,"LA_tree_dist.pickle"), "rb")
        self.node_xy = pickle.load(f1)
        self.adj = pickle.load(f2)
        self.mnodes = pickle.load(f3)
        wyes = pickle.load(f4)
        pred = pickle.load(f5)
        dist = pickle.load(f6)
        f1.close()
        f2.close()
        f3.close()
        f4.close()
        f5.close()
        f6.close()

        is_node_wye = collections.defaultdict(bool)
        for i in wyes:
            is_node_wye[i] = True

        self.nodes = list(self.adj.keys())
        self.weights = {}
        for i in self.nodes:
            if is_node_wye[i]:
                self.weights[i] = 1
            else:
                self.weights[i] = 0

        self.edges = [(i,j) for i in self.nodes for j in self.adj[i] if self.adj[i]!=[]]
        sub_mnodes = list(set(self.nodes).intersection(set(self.mnodes)))
        print(f"Total nodes {len(self.nodes)} total mnodes {len(sub_mnodes)}")
        cut_set, check = self.dynamic_way(self.nodes, self.adj, dist, self.edges)

        self.wye_set = {}
        for i in self.edges:
            self.wye_set[i] = [j for j in cut_set[i] if self.weights[j] > 0]

        t = time.time()
        self.mpots = collections.defaultdict(int)
        for i in sub_mnodes:
            if self.adj[i] != []:
                for j in self.adj[i]:
                    self.mpots[i] += len(self.wye_set[(i, j)])
            else:
                self.mpots[i] = 0

        self.npots = collections.defaultdict(int)
        for i in self.nodes:
            if self.adj[i] != []:
                for j in self.adj[i]:
                    self.npots[i] += len(self.wye_set[i, j])
            else:
                self.npots[i] = 0

        print("starting some unnecessary shallow copying")
        t=time.time()
        self.mnodes=sub_mnodes[:]
        print("time to create shallow copy of variables ",time.time()-t)

    ### send weights code from outside
    @staticmethod
    def truncate_network_using_edges_target_focused(nodes, adj, mnodes, base, edges, npots, cluster_select, weights=None):
        """
        Code edited to match with the fast version
        :param nodes:
        :param adj:
        :param mnodes:
        :param base:
        :param edges:
        :param npots:
        :param cluster_select:
        :param weights:
        :return:
        """
        t = time.time()
        ## npots here is manhole potentials
        if weights==None:
            weights = {}
            for i in nodes:
                if adj[i] == [] and cluster_select[i]:
                    weights[i] = 1
                else:
                    weights[i] = 0
        else:
            for i in nodes:
                if not cluster_select[i]:
                    weights[i]=0
            print(f"Target {sum(weights.values())}")

        to_remove = {i: True for i in edges}

        is_mnode = collections.defaultdict(bool)
        for i in mnodes:
            is_mnode[i] = True

        for i in edges:
            if is_mnode[i[0]] and is_mnode[i[1]]:
                if npots[i[0]] != npots[i[1]] and npots[i[1]] > 0:
                    to_remove[i] = False
            if is_mnode[i[1]] and not is_mnode[i[0]]:
                to_remove[i] = False
            if i[0] == base:
                to_remove[i] = False


        if len(to_remove) > 0:
            edges = [i for i in edges if to_remove[i]]

        # edges=list(set(edges).difference(set(to_remove)))
        #print("Time to get list of edges: ", round(time.time() - t, 20))
        print("Starting iteration", len(edges))
        while len(edges) > 0:
            to_change = edges[0]
            weights[to_change[0]] += weights[to_change[1]]
            adj[to_change[0]].remove(to_change[1])
            # print("To change ",to_change)
            # print("Adj ",adj[to_change[1]])
            for i in adj[to_change[1]]:
                to_remove[(to_change[0], i)] = False
                if to_remove[(to_change[1], i)]:
                    edges.remove((to_change[1], i))
                    # to_remove[(to_change[1], i)]=False
                    if is_mnode[to_change[0]] and is_mnode[i]:
                        if npots[to_change[0]] == npots[i] or npots[i] == 0:
                            edges.append((to_change[0], i))
                            to_remove[(to_change[0], i)] = True

                    else:
                        edges.append((to_change[0], i))
                        to_remove[(to_change[0], i)] = True

                else:
                    if is_mnode[to_change[0]] and is_mnode[i]:
                        if npots[to_change[0]] == npots[i] or npots[i] == 0:
                            edges.append((to_change[0], i))
                            to_remove[(to_change[0], i)] = True

                adj[to_change[0]].append(i)

            nodes.remove(to_change[1])
            if is_mnode[to_change[1]]:
                mnodes.remove(to_change[1])
            adj.pop(to_change[1])
            weights.pop(to_change[1])
            edges.remove(to_change)

        new_edges = [(i, j) for i in nodes for j in adj[i]]
        print(f"Truncation time {time.time() - t}")

        return weights, nodes, adj, mnodes, new_edges

    ## used function
    def get_target_cluster_network(self):
        store_new_files = False

        clusters=self.clusters
        target=self.target[:]
        if self.name=="LA":
            labels, points = self.get_los_angeles_clusters(clusters)
            picdir = la_dir
        elif self.name=="regina":
            labels, points = self.get_regina_clusters(clusters)
            picdir = regina_dir
        cluster_select = collections.defaultdict(bool)
        for i in labels.keys():
            if labels[i] in target:
                cluster_select[i] = True
        self.cluster_select=cluster_select.copy()

        try:
            print(f"Files found for {self.name}")
            f1 = open(os.path.join(picdir,f"{self.name}_trunc_c{target}_adj.pickle"), 'rb')
            f2 = open(os.path.join(picdir,f"{self.name}_trunc_c{target}_mnodes.pickle"), "rb")
            f3 = open(os.path.join(picdir,f"{self.name}_trunc_c{target}_weights.pickle"), "rb")
            f4 = open(os.path.join(picdir,f"{self.name}_trunc_c{target}_mpots.pickle"),"rb")

            self.trunc_adj = pickle.load(f1)
            self.trunc_mnodes = pickle.load(f2)
            self.trunc_weights = pickle.load(f3)
            self.trunc_nodes = list(self.trunc_adj.keys())
            self.tr_mpots=pickle.load(f4)
            self.trunc_edges = [(i,j) for i in self.trunc_nodes for j in self.trunc_adj[i] if self.trunc_adj[i]!=[]]

            f1.close()
            f2.close()
            f3.close()
            f4.close()
            print(f"Target household : {sum([self.trunc_weights[i] for i in self.trunc_weights.keys()])}")
            print(f"MPOTS WS : {self.tr_mpots[self.base]}")
            print(f"Len trunc mnodes : {len(self.trunc_mnodes)}")

            
        except Exception as e:
            print(f"Files not found, doing again and timing for {self.name}. This is the error: ",e)
            
            mcluster = []
            tnodes=set()
            for i in target:
                mcluster.extend([j for j in self.mnodes if labels[j] == i])
                tnodes.update([j for j in labels.keys() if labels[j]==i])

            print(f"Total target nodes {len(tnodes)}")

            
            mpots = collections.defaultdict(int)
            for i in self.nodes:
                if self.adj[i] != []:
                    for j in self.adj[i]:
                        mpots[i] += sum([1 for k in self.wye_set[i, j] if cluster_select[k]])
                else:
                    mpots[i] = 0

            new_w=self.weights.copy()
            for i in self.nodes:
                if not cluster_select[i]:
                    new_w[i]=0

            ## mark important mnodes
            important = []
            for i in self.mnodes:
                if cluster_select[i]:
                    important.append(i)
                else:
                    if self.adj[i] != []:
                        for j in self.adj[i]:
                            if cluster_select[j]:
                                important.append(i)
                                break


            print("Length of nodes and edges before truncation", len(self.nodes), len(self.edges), len(self.adj.keys()))
            print(f"Mnodes before truncation {len(mcluster)}")
            Houses = len([i for i in self.adj.keys() if self.weights[i]>0 and cluster_select[i]])
            print(f"Houses {Houses}, Base {mpots[self.base]}")
            print("Starting truncation")
            t = time.time()
            n=self.nodes[:]
            a=copy.deepcopy(self.adj)
            m=self.mnodes[:]
            e=self.edges[:]
            mp=mpots.copy()

            for i in tnodes:
                plt.plot(self.node_xy[i][0],self.node_xy[i][1],"b.")
            plt.show()

            w=new_w.copy()
            b, nodes, adj, mnodes, edges = ReadData.truncate_network_using_edges_target_focused(n, a, m, self.base, e, mp, cluster_select, w)
            print(f"After trunction nodes {len(nodes)}, edges {len(edges)}, mnodes {len(mnodes)}")
            F = sum([b[i] for i in b.keys()])
            print(max([b[i] for i in mnodes]))
            print("did i cover all houses? ", Houses == F, Houses, F)
            
            if store_new_files:
                f2=open(os.path.join(picdir,f"{self.name}_trunc_c{target}_adj.pickle"),"wb")
                f3=open(os.path.join(picdir,f"{self.name}_trunc_c{target}_mnodes.pickle"),"wb")
                f6=open(os.path.join(picdir,f"{self.name}_trunc_c{target}_weights.pickle"),"wb")
                f7=open(os.path.join(picdir,f"{self.name}_trunc_c{target}_mpots.pickle"),"wb")

                pickle.dump(adj, f2, protocol=pickle.HIGHEST_PROTOCOL)
                pickle.dump(mnodes, f3, protocol=pickle.HIGHEST_PROTOCOL)
                pickle.dump(b, f6, protocol=pickle.HIGHEST_PROTOCOL)
                pickle.dump(mpots,f7,protocol=pickle.HIGHEST_PROTOCOL)
            
                f2.close()
                f3.close()
                f6.close()
                f7.close()

            self.trunc_adj = copy.deepcopy(adj)
            self.trunc_mnodes = mnodes[:]
            self.trunc_weights = b.copy()
            self.trunc_nodes = list(self.trunc_adj.keys())
            self.trunc_edges = [(i,j) for i in self.trunc_nodes for j in self.trunc_adj[i] if self.trunc_adj[i]!=[]]
            self.tr_mpots=mpots.copy()

            print(f"MPOTS WS after new writing: {mpots[self.base]}")


        #self.plot_nodexy(self.trunc_mnodes)

    ## this is for truncating the big network
    def get_trunc_network(self):
        store_new_files = False
        if self.name == 'LA':
            picdir = la_dir
        if self.name == 'regina':
            picdir = regina_dir

        try:
            print(f"Trying to read file for {self.name}")
            f1 = open(os.path.join(picdir,f"{self.name}_trunc_adj.pickle"), 'rb')
            f2 = open(os.path.join(picdir,f"{self.name}_trunc_mnodes.pickle"), "rb")
            f3 = open(os.path.join(picdir,f"{self.name}_trunc_weights.pickle"), "rb")

            self.trunc_adj = pickle.load(f1)
            self.trunc_mnodes = pickle.load(f2)
            self.trunc_weights = pickle.load(f3)
            self.trunc_nodes = list(self.trunc_adj.keys())
            self.trunc_edges = [(i,j) for i in self.trunc_nodes for j in self.trunc_adj[i] if self.trunc_adj[i]!=[]]
            self.tr_mpots=self.mpots.copy()

            f1.close()
            f2.close()
            f3.close()

            print(f"Nodes {len(self.trunc_nodes)}, Manholes {len(self.trunc_mnodes)}, Edges {len(self.trunc_edges)} after truncation")

        except:
            print(f"Files not found, going to write for {self.name}")
            n=self.nodes[:]
            a=copy.deepcopy(self.adj)
            m=self.mnodes[:]
            e=self.edges[:]
            mp=self.mpots.copy()
            w=self.weights.copy()
            TFC_b = sum([w[i] for i in w.keys()])
            w2, n2, a2, m2, e2 = ReadData.truncate_network_using_edges_fast(n, a, m, self.base, e, mp, w)
            print(f"Nodes {len(n2)}, Manholes {len(m2)}, Edges {len(e2)} after truncation")
            TFC = sum([w2[i] for i in w2.keys()])
            print("All houses covered? ", TFC, TFC_b == TFC)

            if store_new_files:
                f2 = open(os.path.join(picdir,f"{self.name}_trunc_adj.pickle"), "wb")
                f3 = open(os.path.join(picdir,f"{self.name}_trunc_mnodes.pickle"), "wb")
                f6 = open(os.path.join(picdir,f"{self.name}_trunc_weights.pickle"), "wb")

                pickle.dump(a2, f2, protocol=pickle.HIGHEST_PROTOCOL)
                pickle.dump(m2, f3, protocol=pickle.HIGHEST_PROTOCOL)
                pickle.dump(w2, f6, protocol=pickle.HIGHEST_PROTOCOL)

                f2.close()
                f3.close()
                f6.close()

            self.trunc_adj = copy.deepcopy(a2)
            self.trunc_mnodes = m2[:]
            self.trunc_weights = w2.copy()
            self.trunc_nodes = list(self.trunc_adj.keys())
            self.trunc_edges = [(i, j) for i in self.trunc_nodes for j in self.trunc_adj[i] if self.trunc_adj[i] != []]
            self.tr_mpots = self.mpots.copy()

            print(f"MPOTS of base after new writing: {self.tr_mpots[self.base]}")

    @staticmethod
    def get_pred_for_mholes(sub_mnodes,adj,base):
        pred={}
        pred[base]=-1
        stack =[base]
        #print("Getting new adjacency for manholes")
        vis=collections.defaultdict(int)
        while len(stack)>0:
            #print(stack)
            curr=stack[0]
            vis[curr]=1
            for i in adj[curr]:
                if vis[i]==0:
                    stack.append(i)
                    if curr==base or curr in sub_mnodes:
                        pred[i]=curr
                    else:
                        pred[i]=pred[curr]
            stack.remove(curr)
        print("Getting adjacency")
        m_adj={}
        for i in sub_mnodes:
            m_adj[i]=[]
            for j in sub_mnodes:
                if i==pred[j]:
                    m_adj[i].append(j)
        m_adj[base]=[i for i in sub_mnodes if pred[i]==base]
        return pred,m_adj

    @staticmethod
    def get_NPS(pots,mnodes,madj,selected):
        nps={}
        for i in mnodes:
            #print("Do BFS to find sensors up stream for selected and non-selected")
            s=[i]
            Ups=[]
            while len(s)>0:
                curr=s[0]
                for k in madj[curr]:
                    if k in selected:
                        Ups.append(k)
                    else:
                        s.append(k)
                s.remove(curr)
            #print(i,Ups)
            if Ups==[]:
                nps[i]=pots[i]
            else:
                nps[i]=pots[i]-sum([pots[k] for k in Ups])
        return nps

    def get_second_layer_reduction(self,S):
        print("Doing second layer reduction directly once")
        F = sum([self.trunc_weights[i] for i in self.trunc_weights.keys()])
        LB = round(F / S)
        pred, madj = ReadData.get_pred_for_mholes(self.trunc_mnodes, self.trunc_adj, self.base)
        nps = ReadData.get_NPS(self.tr_mpots, self.trunc_mnodes, madj, self.trunc_mnodes)
        medges = [(i, j) for i in self.trunc_mnodes for j in madj[i]]
        npred = {}
        for i in self.trunc_nodes:
            for j in self.trunc_adj[i]:
                npred[j] = i
        self.npred=npred.copy()
        G = nx.DiGraph()
        is_mnode = {i: False for i in self.trunc_nodes}
        for i in self.trunc_mnodes:
            is_mnode[i] = True
        n0 = [('s', {"type": "source"})]
        n1 = [(base, {"type": "base", "b": self.trunc_weights[self.base]})]
        n2 = [(i, {"type": "M", "mpots": self.tr_mpots[i], "nps": nps[i], "b": self.trunc_weights[i]})
              for i in self.trunc_mnodes]
        n3 = []
        for i in self.trunc_nodes:
            if i != base and not is_mnode[i]:
                n3.append((i, {"type": "O", "b": self.trunc_weights[i]}))

        G.add_nodes_from(n0)
        G.add_nodes_from(n1)
        G.add_nodes_from(n2)
        G.add_nodes_from(n3)
        e_to_add = [(i, j) for i in self.trunc_nodes for j in self.trunc_adj[i]]
        s_to_m=[('s',i) for i in self.trunc_mnodes]
        s_to_b=[('s',self.base)]
        G.add_edges_from(e_to_add)
        G.add_edges_from(s_to_m)
        G.add_edges_from(s_to_b)
        boundary1 = {}
        boundary2 = {}
        boundary0 = {}

        for (i, j) in medges:
            if self.tr_mpots[i] > LB and self.tr_mpots[j] <= LB:
                boundary0[i] = self.tr_mpots[i]
                boundary1[j] = self.tr_mpots[j]
                for k in madj[j]:
                    boundary2[k] = self.tr_mpots[k]

        for i in madj[base]:
            if self.tr_mpots[i] <= LB:
                boundary1[i] = self.tr_mpots[i]
                for k in madj[i]:
                    boundary2[k] = self.tr_mpots[k]

        temp = list(boundary2.keys())
        red = temp[:]
        while len(temp) > 0:
            e = temp[0]
            for j in madj[e]:
                temp.append(j)
                red.append(j)
            temp.remove(e)

        if remove_red_directly:
            oldG = nx.DiGraph(G)
            for i in boundary1.keys():
                G.nodes[i]['nps'] = G.nodes[i]['mpots']
                G.nodes[i]['b'] = G.nodes[i]['mpots']
            G.remove_nodes_from(red)
            print(f"Number of nodes removed are {len(red)}")
        mnodes = [i for i in G.nodes if G.nodes[i]["type"] == "M"]
        print(f"Number of manholes are {len(mnodes)}")
        adj=nx.to_dict_of_lists(G)
        self.G=nx.DiGraph(G)
        self.b0=copy.deepcopy(boundary0)
        self.b1=copy.deepcopy(boundary1)
        self.b2=copy.deepcopy(boundary2)
        self.red=red[:]
        pred, madj = ReadData.get_pred_for_mholes(mnodes, adj, self.base)
        self.madj=copy.deepcopy(madj)

def FindInBtw(G,thres,madj):
    '''
    gets the inbetween manhole nodes for a valid inequality based on the threshold provided
    '''
    mnodes=[i for i in G.nodes if G.nodes[i]["type"]=="M"]
    minbtw={}
    minbtw_sum={}
    for i in mnodes:
        if G.nodes[i]["mpots"]>thres:
            rem=thres-G.nodes[i]["nps"]
            temp=[i]
            inbtw=[]
            ctr=0
            while len(temp)>0:
                e=temp[0]
                for j in madj[e]:
                    if rem>=G.nodes[j]["nps"]:# and ctr==0:
                        rem=rem-G.nodes[j]["nps"]
                        temp.append(j)
                        inbtw.append(j)
                    elif ctr==0:
                        rem = rem - G.nodes[j]["nps"]
                        temp.append(j)
                        inbtw.append(j)
                        ctr+=1

                temp.remove(e)
            minbtw[i]=inbtw[:]
            minbtw_sum[i]=sum([G.nodes[k]["nps"] for k in inbtw])
        
    return minbtw,minbtw_sum

def primal_model(xvals,nodes,adj):
    t=time.time()
    #print("Solving primal model")
    try:
        inst = cplex.Cplex()
        inst.set_log_stream(None)
        inst.set_error_stream(None)
        inst.set_warning_stream(None)
        inst.set_results_stream(None)
        inst.objective.set_sense(inst.objective.sense.minimize)
        inst.parameters.clocktype.set(2)  # clocktype 2 measures in seconds
        inst.parameters.timelimit.set(6000)  # set the time limit in wall-clock time (in seconds)
        inst.parameters.workmem.set(1024)  # limiting working memory to 1GB, after that node files will be created and stored in HDD
        #inst.parameters.mip.limits.treememory.set(4 * 1024)
        inst.parameters.emphasis.memory.set(1)  # for efficiently using memory
        inst.parameters.workdir.set("NF")  # for storing node files
        inst.parameters.preprocessing.presolve.set(inst.parameters.preprocessing.presolve.values.off)
        #inst.parameters.mip.strategy.file.set(2)  # for storing node files in HDD
        #inst.parameters.mip.pool.capacity.set(1)
        # prob.parameters.mip.strategy.search.set(1)  # 1: B&C, 2: dynamic, 0:let cplex choose
        #inst.parameters.mip.strategy.nodeselect.set(1)  #
        inst.variables.add(obj=[1], lb=[LB], ub=[F], types='C', names=[fmax])
        #for i in f.keys():
        inst.variables.add(obj=[0]*len(f.keys()), lb=[0.0]*len(f.keys()), ub=[F]*len(f.keys()), types='C'*len(f.keys()), names=[f[i] for i in f.keys()])


        inst.linear_constraints.add(lin_expr=[cplex.SparsePair([f['s', i], fmax], [1, -1]) for i in mnodes], senses="L"*len(mnodes), rhs=[0]*len(mnodes), names=[f"c1,{i}" for i in mnodes])
        inst.linear_constraints.add(lin_expr=[cplex.SparsePair([f['s', i]], [1]) for i in mnodes], senses="L"*len(mnodes), rhs=[xvals[i]*mpots[i] for i in mnodes], names=[f"c2,{i}" for i in mnodes])
        inst.linear_constraints.add(lin_expr=[cplex.SparsePair([f['s', i]], [1]) for i in mnodes], senses="G"*len(mnodes), rhs=[xvals[i]*G.nodes[i]["nps"] for i in mnodes], names=[f"c22,{i}" for i in mnodes])


        inst.linear_constraints.add(lin_expr=[cplex.SparsePair([f['s', base]], [1])], senses="E", rhs=[0], names=["c3"])

        for i in mnodes:
            inst.linear_constraints.add(lin_expr=[cplex.SparsePair([f[npred[i], i]], [1])], senses="L", rhs=[F*(1-xvals[i])], names=["c5," + str(i)])
            # if mpots[i] > 0:
            #     inst.linear_constraints.add(lin_expr=[cplex.SparsePair([f[npred[i], i]], [1])], senses="G", rhs=[1-xvals[i]], names=[f"c5p1,{i}"])


        pred = {i: [] for i in trans_nodes}
        for i in adj.keys():
            for j in adj[i]:
                pred[j].append(i)
        for i in nodes:
            if i != 's':
                var = [f[i, j] for j in adj[i]]
                coeff = [-1 for j in adj[i]]
                var.extend([f[j, i] for j in pred[i]])
                coeff.extend([1 for j in pred[i]])
                inst.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="E", rhs=[G.nodes[i]["b"]], names=["c6," + str(i)])
        inst.set_problem_type(inst.problem_type.LP)
        inst.solve()
        #print("primal prob solution status ", inst.solution.get_status_string())
        sol = {i: inst.solution.get_values(i) for i in inst.variables.get_names()}
        obj = inst.solution.get_objective_value()
        gpars.primal_times[0]+=1
        gpars.primal_times[1]+=time.time()-t
        return obj,sol
    except CplexError as e:
        print(f"cplex error {e}")
        return 0,False

def dual_model_docplex(xvals,G):
    print("Dual solver docplex version")
    t=time.time()
    is_mnode={i:True if G.nodes[i]['type']=='M' else False for i in G.nodes}
    mnodes = [i for i in G.nodes if G.nodes[i]['type']=='M']
    try:
    
        inst = dp.Model(name="Dual Model")
        #print("Hi")

        npred = {b:a for (a,b) in G.edges if a!='s'}
        nset = [(j,i) for (i,j) in npred.items()]

        u_base = inst.continuous_var(lb=-cplex.infinity,ub=cplex.infinity,name="u_base")
        u_obj = inst.continuous_var_dict(mnodes,lb=0,ub=cplex.infinity,name="u_obj")
        u_ub = inst.continuous_var_dict(mnodes,lb=0,ub=cplex.infinity,name="u_ub")
        u_lb = inst.continuous_var_dict(mnodes,lb=0,ub=cplex.infinity,name="u_lb")
        u_lbj = inst.continuous_var_dict(nset,lb=0,ub=cplex.infinity,name="u_lb")
        u_ubj = inst.continuous_var_dict(nset,lb=0,ub=cplex.infinity,name="u_ub")
        u_n = inst.continuous_var_dict([i for i in G.nodes if i!='s'],lb=-cplex.infinity,ub=cplex.infinity,name="u_n")
        
        ## dual for fmax
        inst.add_constraint(inst.sum(u_obj)<=1,"fmax_dual")
        ## dual for f_{s,base)
        inst.add_constraint(u_base+u_n[base]<=0,"fsbase_dual")
        ## dual for f_(s,i)
        inst.add_constraints((-u_obj[i]+u_lb[i]-u_ub[i]+u_n[i]<=0 for i in mnodes),"f_s_i")
        ## dual for f_(j,i)
        inst.add_constraints((u_lbj[j,i]-u_ubj[j,i]+u_n[i]-u_n[j]<=0 for (j,i) in nset),"f_j_i")
        ## rest of flow balance dual
        new_set=[(i,j) for i in G.nodes for j in G.neighbors(i) if i!='s' and not is_mnode[j]]
        inst.add_constraints((u_n[j]-u_n[i]<=0 for (i,j) in new_set),"frest_dual")
        ## objective
        objt1= inst.sum(u_lb[i]*xvals[i]*G.nodes[i]["nps"] for i in mnodes)
        objt2= inst.sum(u_ub[i]*xvals[i]*G.nodes[i]["mpots"] for i in mnodes)
        objt3= inst.sum(u_lbj[npred[i],i]*(1-xvals[i])*G.nodes[i]["nps"] for i in mnodes)
        objt4= inst.sum(u_ubj[npred[i],i]*(1-xvals[i])*G.nodes[i]["mpots"] for i in mnodes)
        objt5= inst.sum(u_n[i]*G.nodes[i]["b"] for i in mnodes if i!='s')
        
        objective = objt1-objt2+objt3-objt4+objt5

        inst.add_kpi(objective,"fmax value")
        inst.maximize(objective)
        inst.solve(log_output=False)
        #inst.print_information()
        #inst.report_kpis()
        #print("Hi")

        if "OPTIMAL" in  inst.solve_status.name:
            print("Dual objective matched, getting coefficients for constraint")
            coeffs = {}
            fixc = 0
            for i in xvals.keys():
                Pb = G.nodes[i]["mpots"]
                Qb = G.nodes[i]["nps"]
                coeffs[i] = Qb*(u_lb[i].solution_value-u_lbj[npred[i],i].solution_value)-Pb*(u_ub[i].solution_value-u_ubj[npred[i],i].solution_value)
                fixc += Qb*u_lbj[npred[i],i].solution_value-Pb*u_ubj[npred[i],i].solution_value

            extr = objt5.solution_value +fixc

            return coeffs,extr,objective.solution_value
        
    except CplexError as e:
        print(f"cplex error {e}")
        return 0,False,False

def verificationSubproblem(data,gpars,task, new_cuts,new_lim):
    S=task[2]

    G=nx.DiGraph(data.G)

    ## variables
    fmax = "fmax"
    f = {}
    for i in G.nodes:
        for j in G.neighbors(i):
            f[i, j] = "f" + str(i) + "," + str(j)
    x = {}

    ## update transnodes and transadj and mnodes
    mnodes = [i for i in G.nodes if G.nodes[i]["type"] == "M"]
    for i in mnodes:
        x[i] = "x" + str(i)

    trans_nodes = [i for i in G.nodes]
    trans_adj = nx.to_dict_of_lists(G)

    try:
        prob=cplex.Cplex()
        prob.objective.set_sense(prob.objective.sense.minimize)
        #prob.parameters.mip.display.set(4)
        prob.set_log_stream(None)
        prob.set_error_stream(None)
        prob.set_warning_stream(None)
        prob.set_results_stream(None)
        prob.parameters.clocktype.set(2)  # clocktype 2 measures in seconds
        prob.parameters.timelimit.set(100)  # set the time limit in wall-clock time (in seconds)
        prob.parameters.workmem.set(4*1024)  # limiting working memory to 1GB, after that node files will be created and stored in HDD
        prob.parameters.mip.limits.treememory.set(10*1024)
        prob.parameters.emphasis.memory.set(1)  # for efficiently using memory
        prob.parameters.workdir.set(nodefile_dir)  # for storing node files
        prob.parameters.mip.strategy.file.set(2)  # for storing node files in HDD
        prob.parameters.preprocessing.presolve.set(0)
        prob.parameters.mip.strategy.nodeselect.set(1)  # 0: DFS, 1: Best Bound S, 2: Best estimate S, 3: Alternative BES, default is 1

        file2 = os.path.join(lpdir,f"{gpars.id}.lp")
        prob.read(file2)
        print("Read file for solving")


        for i in new_cuts.keys():
            expr=[cplex.SparsePair([x[j] for j in new_cuts[i]],[1 for j in new_cuts[i]])]
            rhs=[1]
            prob.linear_constraints.add(lin_expr=expr,senses="G",rhs=rhs,names=[f"new_cuts_{i}"])

        print("Added cuts")
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(["fmax"],[1])],senses="L",rhs=[new_lim],names=["obj_limit"])
        print("Solving for feasibility")
        prob.parameters.mip.limits.solutions.set(1)

        if use_new_vi:
            aa,bb=FindInBtw(G,new_lim,data.madj)
            print(f"Adding additional VI for {len(aa)} manholes inbtw")
            for i in aa.keys():
                expr=[cplex.SparsePair([x[j] for j in aa[i]],[1 for j in aa[i]])]
                rhs=[1]
                prob.linear_constraints.add(lin_expr=expr,senses="G",rhs=rhs,names=[f"new_vi_cuts_{i}"])



        t=time.time()
        prob.solve()
        t2=time.time()
        print(f"Status : {prob.solution.get_status_string()}")

        infeas=[103,106,108,112]
        if prob.solution.get_status() not in infeas:
            flow={}
            for i in f.keys():
                flow[i]=round(prob.solution.get_values(f[i]),3)

            xv={}
            for i in mnodes:
                #print("Selected ",x[i],prob.solution.get_values(x[i]))
                xv[i]=round(prob.solution.get_values(x[i]),3)
            tf=0

            selected=[]
            flt=[]
            fm=0
            fs=0
            for i in mnodes:
                if xv[i]>0:
                    if flow['s',i]>fm:
                        fm=flow['s',i]
                        fs=i
                    selected.append(i)
                    flt.append((i,flow['s',i]))
                    #print("Flow {0} {1} Selected {2}".format(flow['s',i],('s',i),xv[i]))
                    tf+=flow['s',i]

            #gpars.save_csv_solution(prob,data,task,selected,flt,Houses,t2-t,flow)
            print("Objective ",prob.solution.get_objective_value())
            print("Total flow covered : ",tf)
            print(f"Sensor with max flow {fs}")
            print("Max flow : ",fm)
            print("Min flow : ",prob.solution.get_values("fmin"))
            return True, fm,selected

        elif prob.solution.get_status()==103:
            print("Infeasible")
            return False,False,0
        else:
            return False,False,1


    except TypeError as e:
        print(f"Errrrr {e}")
        return False,False,2

## this is being used currently
def flowCoveringLocationProblem(data,gpars,task, thres=False, new_lb=0, best_inc=0):
    '''
    task = [strategy, model, number of sensors, no of clusters (optional), target cluster (optional)]
    '''
    S=task[2]
    strategy=task[0]

    print(f"Solving strategy {strategy}")
    nodes=data.trunc_nodes[:]
    adj=copy.deepcopy(data.trunc_adj)
    mnodes=data.trunc_mnodes[:]
    base=data.base
    edges=data.trunc_edges[:]
    mpots=data.tr_mpots.copy()
    npots=data.npots.copy()
    weights=data.trunc_weights.copy()
    print("Length of nodes and edges and mnodes ",len(nodes),len(edges),len(adj.keys()),len(mnodes))
    if weights==None:
        Houses = sum([weights[i] for i in weights.keys()])
    else:
        Houses=sum([weights[i] for i in weights.keys()])
    print("Starting IP")
    t=time.time()
    b= weights.copy()
    LB=Houses/S
    
    ### quick BFS check to see if all nodes are reachable after truncation - optional
    print("Quick BFS Check")
    t=time.time()
    is_visited={i:False for i in nodes}
    temp=[base]
    is_visited[base]=True
    while len(temp)>0:
        temp1=temp[0]
        for i in adj[temp1]:
            if not is_visited[i]:
                temp.append(i)
                is_visited[i]=True
        temp.remove(temp1)


    npred=data.npred.copy()
    #print("Length of nodes and edges", len(nodes), len(edges), len(adj.keys()),time.time()-t)
    F=sum([b[i] for i in b.keys()])
    #print(max([b[i] for i in mnodes]))
    print("did i cover all houses? ",Houses==F,Houses,F)

    ### create a directed graph to better handle stuff
    G = nx.DiGraph(data.G)
    mnodes = [i for i in G.nodes if G.nodes[i]["type"] == "M"]
    is_mnode={i:False for i in G.nodes}
    for i in mnodes:
        is_mnode[i]=True

    fmax = "fmax"
    f = {}
    for i in G.nodes:
        for j in G.neighbors(i):
            f[i, j] = "f" + str(i) + "," + str(j)
    x = {}

    for i in mnodes:
        x[i] = "x" + str(i)

    trans_nodes = [i for i in G.nodes]
    trans_adj = nx.to_dict_of_lists(G)

    print(f"After direct removal of redundant nodes {second_layer_red}: {remove_red_directly}")
    print("Length of nodes in transformed network ", len(trans_nodes), len(trans_adj.keys()))
    print(f"Number of mnodes : {len(mnodes)}, base is : {base}")
    #pred, madj = ReadData.get_pred_for_mholes(mnodes, trans_adj, base)
    madj=copy.deepcopy(data.madj)
    # data.store_network_used(madj,boundary0,boundary1,boundary2)
    ### implements lazy cuts new using context

    class GenericCB():
        is_edge_used=collections.defaultdict(bool)
        is_node_used=collections.defaultdict(bool)
        use_cuts=False
        use_new_vi_cuts=False
        new_bound=LB
        use_new_bound=False
        cut_pool={}
        vi_cut_pool={}
        relax_invoked=0
        couldnt=set([Houses])
        def invoke(self, context):
            if context.in_candidate():# or context.in_global_progress():
                flow = {}
                for i in f.keys():
                    flow[i] = round(context.get_candidate_point(f[i]), 3)

                xv = {}
                selected = []
                flt = []
                for i in mnodes:
                    xv[i] = round(context.get_candidate_point(x[i]), 3)
                    if xv[i] > 0:
                        selected.append(i)
                        flt.append((i, context.get_candidate_point(f['s', i])))

                fmax_val = round(context.get_candidate_point("fmax"))

                best_inc_bef = round(context.get_incumbent_objective())

                curr_inc = round(context.get_candidate_objective())
                print(f"Invoked in candidate context {curr_inc}, relax invoked {self.relax_invoked}")

                best_bound = context.get_double_info(context.info.best_bound)

                best_inc=curr_inc

                case = [(j, i) for (j, i) in G.edges if G.nodes[i]["type"] == "M" and G.nodes[j]["type"] == "M" and
                         mpots[j] > best_inc and mpots[i] <= best_inc and mpots[i] > 0 and not self.is_edge_used[j,i,best_inc]]

                treej = {}
                treesumj = {}
                flag=False
                for (j, i) in case:
                    self.is_edge_used[j,i,best_inc]=True
                    if not self.is_node_used[j,best_inc]:
                        self.is_node_used[j,best_inc]=True
                        treej[j] = []
                        temp = [j]
                        while len(temp) > 0:
                            curr = temp[0]
                            for k in madj[curr]:
                                if mpots[k] <=best_inc and mpots[k] > 0:
                                    treej[j].append(k)
                                temp.append(k)
                            temp.remove(curr)
                        treesumj[j] = sum([xv[k] for k in treej[j]])
                        self.cut_pool[j,best_inc]=[treej[j][:],True]
                        flag=True
                        self.use_cuts=True
                #print("Hi")
                if flag:
                    print(f"Some cuts queued for obj {best_inc}")

                if use_new_vi:
                    aa,bb=FindInBtw(G,best_inc,madj)
                    self.use_new_vi_cuts=True
                    self.vi_cut_pool=copy.deepcopy(aa)

            ## is not being used
            if context.in_relaxation() and self.use_new_bound:
                best_bound = context.get_double_info(context.info.best_bound)
                if best_bound<self.new_bound:
                    self.relax_invoked += 1
                    #print(f"Imposing new bound {self.new_bound}")
                    context.add_user_cuts(cuts=[cplex.SparsePair(["fmax"],[1])],senses={"G"},rhs=[self.new_bound],cutmanagement=[0],local=[0])

            if context.in_relaxation() and self.use_new_vi_cuts:
                print("Adding new VI cuts")
                self.relax_invoked+=1
                for i in self.vi_cut_pool.keys():
                    expr=[cplex.SparsePair([x[j] for j in self.vi_cut_pool[i]],[1 for j in self.vi_cut_pool[i]])]
                    rhs=[1]
                    context.add_user_cuts(cuts=expr,senses="G",rhs=rhs,cutmanagement=[0],local=[0])
                self.use_new_vi_cuts=False

            if context.in_relaxation() and self.use_cuts:
                if len(self.cut_pool.keys())>0:
                    print("Adding cuts based on current inc and and anything below")
                    self.relax_invoked+=1
                    for i in self.cut_pool.keys():
                        if self.cut_pool[i][1]:
                            expr=cplex.SparsePair([x[j] for j in self.cut_pool[i][0]],[1 for j in self.cut_pool[i][0]])
                            rhs=1
                            context.add_user_cuts(cuts=[expr],senses={"G"},rhs=[rhs],cutmanagement=[0],local=[0])
                            self.cut_pool[i][1]=False
                    self.cut_pool={}

                #print("Hi")

    try:
        prob=cplex.Cplex()
        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.parameters.clocktype.set(2)  # clocktype 2 measures in seconds
        prob.parameters.timelimit.set(gpars.time_limit)  # set the time limit in wall-clock time (in seconds)
        prob.parameters.workmem.set(4*1024)  # limiting working memory to 1GB, after that node files will be created and stored in HDD
        prob.parameters.mip.limits.treememory.set(10*1024)
        prob.parameters.emphasis.memory.set(1)  # for efficiently using memory
        prob.parameters.workdir.set(nodefile_dir)  # for storing node files
        prob.parameters.mip.strategy.file.set(2)  # for storing node files in HDD
        prob.parameters.mip.tolerances.absmipgap.set(0.99999)
        prob.parameters.mip.strategy.nodeselect.set(1)  # 0: DFS, 1: Best Bound S, 2: Best estimate S, 3: Alternative BES, default is 1
        fmin="fmin"
        if "lb" in strategy:
            prob.variables.add(obj=[1], lb=[new_lb], ub=[F], types='C', names=[fmax])
        else:
            prob.variables.add(obj=[1], lb=[LB], ub=[F], types='C', names=[fmax])
        prob.variables.add(obj=[0],lb=[0],ub=[LB],types='C',names=[fmin])
        
        for i in f.keys():
            if i[0]=='s':
                prob.variables.add(obj=[0], lb=[0.0], ub=[F], types='C', names=[f[i]])
            else:
                prob.variables.add(obj=[0], lb=[0.0], ub=[F], types='C', names=[f[i]])
        for i in x.keys():
            prob.variables.add(obj=[0], lb=[0], ub=[1], types='I', names=[x[i]])
        print("Adding constraints")
        
        ## 1
        print("Adding constraint 1")
        for i in mnodes:
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair([f['s',i],fmax],[1,-1])],senses="L",rhs=[0],names=["fmax,"+str(i)])
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair([f['s',i], x[i]], [1, -mpots[i]])], senses="L", rhs=[0],names=["c2," + str(i)])
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair([f['s',i], x[i]], [1, -G.nodes[i]["nps"]])],senses="G", rhs=[0], names=["c22," + str(i)])
            
        ## 3
        print("Adding constraint 3")
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair([f['s', base]], [1])], senses="E", rhs=[0], names=["c3"])
        
        ## 4
        print("Adding constraint 4")
        var = []
        coeff = []
        for i in mnodes:
            var.append(x[i])
            coeff.append(1)
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="E", rhs=[S], names=["c4"])


        ## 5
        print("Adding constraint 5")
        for i in mnodes:
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair([f[npred[i], i],x[i]], [1,mpots[i]])], senses="L", rhs=[mpots[i]], names=["c5,"+str(i)])
            if mpots[i]>0:
                prob.linear_constraints.add(lin_expr=[cplex.SparsePair([f[npred[i], i], x[i]], [1, G.nodes[i]["nps"]])], senses="G",rhs=[G.nodes[i]["nps"]], names=[f"c5p1,{i}"])
            else:
                prob.linear_constraints.add(lin_expr=[cplex.SparsePair([x[i]], [1])], senses="E",rhs=[0], names=[f"c5p2,{i}"])

        ## 6
        print("Adding constraint 6")
        trans_pred={i:[] for i in trans_nodes}
        for i in trans_adj.keys():
            for j in trans_adj[i]:
                trans_pred[j].append(i)
        for i in trans_nodes:
            if i!='s':
                var=[f[i,j] for j in trans_adj[i]]
                coeff=[-1 for j in trans_adj[i]]
                var.extend([f[j,i] for j in trans_pred[i]])
                coeff.extend([1 for j in trans_pred[i]])
                prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="E", rhs=[G.nodes[i]["b"]], names=["c6,"+str(i)])
            else:
                var = [f[i, j] for j in trans_adj[i]]
                coeff = [-1 for j in trans_adj[i]]
                var.extend([f[j, i] for j in trans_pred[i]])
                coeff.extend([1 for j in trans_pred[i]])
                prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="E", rhs=[-F], names=["c6," + str(i)])

        print("Starting to solve")

        def add_flow_vis(mnodes,madj,mpots):
            for i in mnodes:
                var=[f["s",str(i)+"d"]]
                coeff=[1]
                s=[i]
                while len(s):
                    curr=s[0]
                    for j in madj[curr]:
                        s.append(j)
                        var.append(f["s",str(j)+"d"])
                        coeff.append(1)
                    s.remove(curr)
                #var.append(x[i])
                #coeff.append(-mpots[i])
                print(var,coeff)
                prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var,coeff)],senses="L",rhs=[mpots[i]],names=["VI1_"+str(i)])


        ## valid inequality idea
        if "basic" not in strategy:
            print("Also adding valid inequalities") ## constraints 4 from paper
            expr=[]
            rhs=[]
            names=[]
            for (j,i) in G.edges:
                if G.nodes[j]["type"]=="M" and G.nodes[i]["type"]=="M" and G.nodes[npred[j]]["type"]=="M":
                    var=[f['s',j],f[npred[j],j],x[i]]
                    coeff=[1,1,G.nodes[i]["nps"]]
                    expr.append(cplex.SparsePair(var,coeff))
                    rhs.append(G.nodes[j]["nps"]+G.nodes[i]["nps"])
                    names.append(f"new_vi_{(j,i)}")
            sens="G"*len(expr)
            prob.linear_constraints.add(lin_expr=expr,senses=sens,rhs=rhs,names=names)

            ### VI2 : set of mnodes whose potential is just below ideal threshold and their pred is 'ws', must always be selected to achieve full coverage
            must_select=[i for i in G.neighbors(base) if mpots[i]<LB and G.out_degree(i)==0 and mpots[i]>0]
            print(f"Must select ",must_select)
            for i in must_select:
                prob.linear_constraints.add(lin_expr=[cplex.SparsePair([x[i]],[1])],senses="E",rhs=[1],names=[f"must_select_vi_{i}"])


        if "lb" in strategy and best_inc==0:
            '''
            this is only executed during contextwithlb
            '''
            thres=np.ceil(thres)
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair(["fmax"],[1])],senses="L",rhs=[thres])
            aa,bb=FindInBtw(G,thres,madj)
            print(f"Adding additional VI for {len(aa)} manholes inbtw and thres of {thres}")
            for i in aa.keys():
                expr=[cplex.SparsePair([x[j] for j in aa[i]],[1 for j in aa[i]])]
                rhs=[1]
                prob.linear_constraints.add(lin_expr=expr,senses="G",rhs=rhs,names=[f"new_vi_cuts_{i}"])

            case = [(j, i) for (j, i) in G.edges if G.nodes[i]["type"] == "M" and G.nodes[j]["type"] == "M" and
                    mpots[j] > thres and mpots[i] <= thres and mpots[i] > 0]

            print(f"Some more constraint {len(case)}")
            treej = {}
            treesumj = {}
            for (j, i) in case:
                treej[j] = []
                temp = [j]
                while len(temp) > 0:
                    curr = temp[0]
                    for k in madj[curr]:
                        if mpots[k] <= thres and mpots[k] > 0:
                            treej[j].append(k)
                        temp.append(k)
                    temp.remove(curr)

                expr = [cplex.SparsePair([x[k] for k in treej[j]], [1 for k in treej[j]])]
                rhs = [1]
                prob.linear_constraints.add(lin_expr=expr, senses="G", rhs=rhs, names=[f"new_cuts_{j}"])
            prob.parameters.mip.limits.solutions.set(1)
            prob.solve()
            infeas = [103, 106, 108, 112]
            if prob.solution.get_status() not in infeas:
                selected=[i for i in mnodes if prob.solution.get_values(x[i])>0]
                return "feas",prob.solution.get_objective_value(),selected
            elif prob.solution.get_status()==103:
                return "infeas",0,[]
            else:
                return "stop",0,[]


        if best_inc>0:
            '''
            this is only executed during contextwithlb
            '''
            best_inc=round(best_inc)
            thres=np.ceil(thres)
            print(f"Adding empirical LB of {thres}")
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair([fmax], [1])], rhs=[thres], senses="G", names=[f"th_lb"])
            aa, bb = FindInBtw(G, best_inc, madj)
            print(f"Adding additional VI for {len(aa)} manholes inbtw and best inc of {best_inc}")
            for i in aa.keys():
                expr = [cplex.SparsePair([x[j] for j in aa[i]], [1 for j in aa[i]])]
                rhs = [1]
                prob.linear_constraints.add(lin_expr=expr, senses="G", rhs=rhs, names=[f"new_vi_cuts_{i}"])

            case = [(j, i) for (j, i) in G.edges if G.nodes[i]["type"] == "M" and G.nodes[j]["type"] == "M" and
                    mpots[j] > best_inc and mpots[i] <= best_inc and mpots[i] > 0]

            print(f"Some more constraint {len(case)}")
            treej = {}
            treesumj = {}
            for (j, i) in case:
                treej[j] = []
                temp = [j]
                while len(temp) > 0:
                    curr = temp[0]
                    for k in madj[curr]:
                        if mpots[k] <= best_inc and mpots[k] > 0:
                            treej[j].append(k)
                        temp.append(k)
                    temp.remove(curr)

                expr = [cplex.SparsePair([x[k] for k in treej[j]], [1 for k in treej[j]])]
                rhs = [1]
                prob.linear_constraints.add(lin_expr=expr, senses="G", rhs=rhs, names=[f"new_cuts_{j}"])
            print(f"Also adding current best inc as warm start {best_inc}")
            varr=[x[i] for i in data.selected]
            vall=[1 for i in data.selected]
            prob.MIP_starts.add([varr, vall], prob.MIP_starts.effort_level.repair, "first")

        if "context" in strategy:
            cm = Context.id.candidate | Context.id.relaxation
            prob.set_callback(GenericCB(), cm)


        file2=os.path.join(lpdir,f"{gpars.id}.lp")
        prob.write(file2)
        t=time.time()
        prob.solve()
        t2=time.time()

        infeas=[103,106,108,112]
        if prob.solution.get_status() not in infeas:
            flow={}
            for i in f.keys():
                flow[i]=round(prob.solution.get_values(f[i]),3)

            xv={}
            for i in mnodes:
                xv[i]=round(prob.solution.get_values(x[i]),3)
            tf=0

            selected=[]
            flt=[]
            for i in mnodes:
                if xv[i]>0:
                    selected.append(i)
                    flt.append((i,flow['s',i]))
                    print("Flow {0} {1} Selected {2}".format(flow['s',i],('s',i),xv[i]))
                    tf+=flow['s',i]
            
            gpars.save_csv_solution(prob,data,task,selected,flt,Houses,t2-t,flow)
            print("Objective ",prob.solution.get_objective_value())
            print("Gap ",prob.solution.MIP.get_mip_relative_gap())
            print("Total flow covered : ",tf)
            print("Max flow : ",prob.solution.get_values(fmax))
            print("Min flow : ",prob.solution.get_values(fmin))
            return flow, xv

        else:
            print("Infeasible")
            gpars.save_csv_solution(prob, data, task, [], [], Houses, t2 - t, [])
            print("Objective ", gpars.current_best)
            return True, True


    except CplexError as e:
        print(f"Errrrr {e}")
        return False,False

def binary_search_algorithm(data,gpars,task,wt):
    print("Running binary search")
    LB = sum(data.trunc_weights.values()) / task[2]
    new_LB = LB
    new_UB = LB * task[2]
    t = time.time()
    l = LB
    u = LB * task[2]


    while u - l > 10:
        thres = (u + wt*l) / (1+wt)
        print(f"trying u {u} l {l} thres {thres}")
        st, obj, sel = flowCoveringLocationProblem(data, gpars, task, thres, l)
        print(f"Thres {thres} Status {st}")
        if st == "stop":
            u = thres

        elif st == "feas":
            u = obj
            new_UB = obj
            data.store_solution(sel)
        else:
            l = thres
    t2=time.time()
    print(f"Total time to find best lower bound {t2 - t}")
    print(f"u {u}, l {l}, best inc {new_UB}")
    new_LB = l
    gpars.save_ws_time(t2-t,new_LB,new_UB)
    return new_LB,new_UB,data,gpars



if __name__=="__main__":
    second_layer_red=True
    use_new_vi=True
    remove_red_directly=True

    print("In reference to section 4.2 of the paper, the four instances are defined as follows.")
    '''
    M1: problem solved without network pre-processing and without any acceleration techniques.
    M2: problem solved with network pre-processing but without any accelation techniques.
    M3: problem solved with network pre-processing and boundary layer reduction.
    M4: problem solved with all reduction techniques and acceleration techniques.
    
    Task description is as follows: (4th and 5th index is only provided for "flow_clus" model).
    orig_task[0]: strategy - M1, M2, M3, M4
    orig_task[1]: model - "flow", "flow_clus"
    orig_task[2]: S - number of sensors (any number more than 10, below that may make the problem infeasible)
    orig_task[3]: target zone for "flow_clus" model. For LA, it can be either of 0,1,2. For regina, it can be either of 1,2,3,4,5.
            target zone must be passed as a list. For example, to run for Zone 3 of regina, orig_task[3] must be [3]

    For example, inst_name = "regina" and orig_task = ["M4","flow_clus",20,[2]] will run the M4 strategy for flow_clus model with 20 sensors for regina with target Zone 2.
    '''

    
    inst_name = "regina"
    orig_task = ["M4","flow_clus",20,[1]]
    save_solution_map = True

    task_generator = create_task(inst_name,orig_task)
    task = task_generator.get_task()
    
    print(f"Processing task {task}, {inst_name}")

    if inst_name == "LA":
        base = 'ws'
    elif inst_name == "regina":
        base = 149390
    else:
        base = 0


    if "lb" in task[0]:
        if "clus" in task[1]:
            gpars = global_params(f"{inst_name}_{task[0]}_{task[1]}_{task[2]}_{task[4]}", 100)
            data = ReadData(name=inst_name, base=base, noc=task[3], target=task[4])
            data.get_target_cluster_network()
        else:
            gpars = global_params(f"{inst_name}_{task[0]}_{task[1]}_{task[2]}", 100)
            data = ReadData(name=inst_name, base=base) 
            data.get_trunc_network()

        data.get_second_layer_reduction(task[2])

        new_LB, new_UB, data, gpars = binary_search_algorithm(data, gpars, task, 1)
        wst = gpars.wst
        if "clus" in task[1]:
            gpars = global_params(f"{inst_name}_{task[0]}_{task[1]}_{task[2]}_{task[4]}", 7200 - gpars.wst)
        else:
            gpars = global_params(f"{inst_name}_{task[0]}_{task[1]}_{task[2]}", 7200 - gpars.wst)
        gpars.save_ws_time(wst, new_LB, new_UB)
        f, x = flowCoveringLocationProblem(data, gpars, task, thres=new_LB, new_lb=new_LB, best_inc=new_UB)
        

    else:
        if "clus" in task[1]:
            gpars = global_params(f"{inst_name}_{task[0]}_{task[1]}_{task[2]}_{task[4]}")
            data = ReadData(name=inst_name, base=base, noc=task[3], target=task[4])
            data.get_target_cluster_network()
        else:
            gpars = global_params(f"{inst_name}_{task[0]}_{task[1]}_{task[2]}")
            data = ReadData(name=inst_name, base=base)  # , noc=task[3], target=task[4])
            data.get_trunc_network()

        if "nored" in task[0]:
            remove_red_directly = False
        data.get_second_layer_reduction(task[2])

        f, x = flowCoveringLocationProblem(data, gpars, task)


    if save_solution_map:
        print("Going to save solution map in results folder")
        data.generate_solution_shapefile(task)

