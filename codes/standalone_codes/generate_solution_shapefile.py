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

folder="./ResultsLatest/"
picdir = "PickleFiles/"


inst="LA"
strat="contextwithlb"
model="flow_clus"
sens=100
target=[0]

filename=folder+f"{inst}_{strat}_{model}_{sens}_{target}.csv"
data=pd.read_csv(filename)

locs=data.Locs.values[0].replace("'"," ")
locs=locs.replace("["," ")
locs=locs.replace("]"," ")
locs=locs.replace(","," ")
LOCS=[int(i) for i in locs.split(" ") if len(i)>0]

f1 = open(picdir + "LA_new_cord.pickle", 'rb')
node_xy = pickle.load(f1)
f1.close()

PL=[Point(node_xy[i]) for i in LOCS]
newdata = geopandas.GeoDataFrame(geometry=PL)
newdata.to_file(folder+f"{inst}_{strat}_{model}_{sens}_{target}.shp")