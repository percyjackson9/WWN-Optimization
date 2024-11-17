#!/usr/bin/env python
import numpy as np

def run0():
    nodes = ["flow","disj"]
    p = [20,50,100]
    task_list=[[a,b] for a in nodes for b in p]
    tasks={}
    for (i,j) in enumerate(task_list):
        tasks[i+1]=j
    return tasks

if __name__ == '__main__':
    import pickle,sys
    pickle_file = "records/"
    with open(pickle_file+"run1_sewage.pykl",'wb') as fh:
        tasks = run0()
        pickle.dump(tasks,fh)
