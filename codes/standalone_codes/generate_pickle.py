import pickle
import pandas as pd

strategy=["oneshot","lazyrej","benders","basic"]
strategy=["context"]#,"contextrej"]#,"benders","basic"]
strategy=["basicnored","basic","oneshot","contextwithlb"]
model=["flow","flow_clus"]
S=[10]
CS=7
clusters=[[0],[1],[2]]

ResultsF="./ResultsNewTry/"

task_list={}
ctr=1
i_list=[]
k_list=[]
j_list=[]
l_list=[]

ii_list=[]
kk_list=[]
jj_list=[]
d1=pd.DataFrame()
d2=pd.DataFrame()

for k in S:
    for j in model:
        if "clus" in j:
            for l in clusters:
                for i in strategy:
                    task_list[ctr]=[i,j,k,CS,l]
                    # fname=ResultsF+f"los_angeles_{i}_{j}_{k}_{l}.csv"
                    # data=pd.read_csv(fname)
                    # data=data.drop(columns=['Unnamed: 0'])
                    # data["Strategy"]=[i]
                    # data["NOS"]=[k]
                    # data["Target"]=[l]
                    # d1=d1.append(data)
                    ctr+=1
        else:
            for i in strategy:
                task_list[ctr]=[i,j,k]
                # fname = ResultsF + f"los_angeles_{i}_{j}_{k}.csv"
                # data = pd.read_csv(fname)
                # data = data.drop(columns=['Unnamed: 0'])
                # data["Strategy"] = [i]
                # data["NOS"] = [k]
                # d2 = d2.append(data)
                ctr+=1
#
# with open("wwn_contextonly_instances.pykl","wb") as f:
#     pickle.dump(task_list,f)
#
# with open("wwn_contextonly_instances.pykl","rb") as f:
#     tasks=pickle.load(f)

# d1.to_csv("wwn_results_flow_clus_context_new.csv")
# d2.to_csv("wwn_results_flow_context_new.csv")

with open("LA_10_sensors_only.pykl","wb") as f:
    pickle.dump(task_list,f)