import pandas as pd
import numpy as np
import os
from pathlib import PurePath
import pickle

curd = PurePath(os.path.dirname(os.path.abspath(__file__)))
pdir = curd.parent
picdir = os.path.join(pdir,"pickle_files")
la_dir = os.path.join(picdir,"la_network")
regina_dir = os.path.join(picdir,"regina_network")
results_dir = os.path.join(pdir,"results")
lp_file_dir = os.path.join(pdir,"lp_files")
lpdir = lp_file_dir
nodefile_dir = lp_file_dir

file_to_read = "LA_trunc_c[1]_mpots.pickle"
with open(os.path.join(la_dir,file_to_read),"rb") as f:
    la = pickle.load(f)