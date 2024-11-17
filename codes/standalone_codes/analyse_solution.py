'''
This is a standalone code not used anywhere else
'''

import numpy as np

def analyse_solution(nodes,adj,selected,base,weights,m_adj):
    edges=[]
    for i in nodes:
        for j in adj[i]:
            edges.append((i,j))

    is_boundary={}
    for i in selected:
        is_boundary[i]=True
        a=[i]
        while len(a)>0:
            s=a[0]
            flag=True
            for j in m_adj[s]:
                a.append(j)
                if j in selected:
                    is_boundary[i]=False
                    flag=False
                    break
            if not flag:
                break
            a.remove(s)
    boundary=[i for i in selected if is_boundary[i]==True]
    cases=np.eye(len(boundary),dtype=int)
    total_unc={}
    for ii,i in enumerate(cases):
        wedges={}
        for j in edges:
            wedges[j] = 0
        no_edges = []
        for j,k in enumerate(boundary):
            if i[j]==0:
                for s in adj[k]:
                    no_edges.append((k,s))

        for j in edges:
            if j in no_edges:
                wedges[j] = 1
        a = []  ## set of reachable nodes
        b = [base]
        while len(b) > 0:
            e = b[0]
            a.append(e)
            if adj[e] != []:
                for j in adj[e]:
                    if wedges[e, j] == 0:
                        b.append(j)
                        a.append(j)
            b.remove(e)

        total_unc[ii] = sum(weights[i] for i in a)
        print("Case,Unc : {0} , {1} ".format(i, total_unc[ii]))
    return total_unc