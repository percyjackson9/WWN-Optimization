'''
Used to visualize solutions in GIS.
'''

from shapely.geometry import Point
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import collections
import time
import pickle
import geopandas
import os
from pathlib import PurePath


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



inst_name="regina"
orig_task = ["M4","flow_clus",20,[1]]
task_generator = create_task(inst_name,orig_task)
task = task_generator.get_task()


filename=f"{inst_name}_{task[0]}_{task[1]}_{task[2]}_{task[4]}.csv"
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
newdata.to_file(os.path.join(results_dir,f"{inst_name}_{task[0]}_{task[1]}_{task[2]}_{task[4]}.shp"))