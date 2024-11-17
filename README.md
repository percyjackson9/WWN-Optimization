Instructions for running the codes for the paper titles "An Optimization Approach for Biosurveillance in Wastewater Networks". The code and datasets are licensed under MIT License. Please cite our article if you use the code and instances!

**Folder Structure**
1. /codes/: this contains the main code named "flow_covering_model_code.py".
2. /lp_files/: this folder stores temporary lp files of the model. The files are safe (and recommended) to delete after a successful run is complete.
3. /pickle_files/la_network/: this folder stores all the relevant data for LA network.
4. /pickle_files/regina_network/: this folder stores all the relevant data for Regina network.
5. /results/: this folder stores the results generated for each run in the form of CSV, some pickle files, and shapefiles.

**LA Network File description**
1. LA_new_cord.pickle: a dictionary for coordinate of the LA network.
2. LA_new_mnodes.pickle: a set of nodes which are considered manholes.
3. LA_new_wyes.pickle: a set of nodes which are considered wyes.
4. LA_tree_adj.pickle: adjacency list of the directed tree network, where each edge direction is opposite to the flow of wastewater.
5. LA_tree_dist.pickle: dictionary, where each key is a node in the directed tree, and value is the distance in 'number' of edges from 'ws' to that node.
6. LA_tree_pred.pickle: oppositive of the adjacency file. Since its a tree, for every key in this dictionary there is a single value representing its parent node.
7. LA_trunc_mnodes.pickle: set of manholes after network pre-processing as described in the paper.
8. LA_trunc_weights.pickle: a dictionary where the key is a manhole node and value is the minimum potential of that manhole (refer to paper).
9. los_angeles_7.pickle: a dictionary that labels each node as a number between 0 and 6 (inclusive), assigning a unique cluster. These clusters are used as Zones.
10. LA_trunc_c[i]_adj.pickle: adjacency list of cluster i after pre-processing.
11. LA_trunc_c[i]_mnodes.pickle: list of manholes in cluster i after pre-processing.
12. LA_trunc_c[i]_mpots.pickle: max potential of each manhole after pre-processing.
13. LA_trunc_c[i]_weights.pickle: min potential of each manhole after pre-processing.


**Regina Network File description**
1. regina_all_cord.pickle: a dictionary for coordinate of the regina network.
2. regina_all_manholes.pickle: a set of nodes which are considered manholes.
3. regina_all_adj.pickle: adjacency list of the directed tree network, where each edge direction is opposite to the flow of wastewater.
4. regina_trunc_adj.pickle: adjacency list of the directed tree network, after preprocessing.
7. regina_trunc_mnodes.pickle: set of manholes after network pre-processing as described in the paper.
8. regina_trunc_weights.pickle: a dictionary where the key is a manhole node and value is the minimum potential of that manhole (refer to paper).
9. regina_10.pickle: a dictionary that labels each node as a number between 0 and 9 (inclusive), assigning a unique cluster. These clusters are used as Zones.
10. regina_trunc_c[i]_adj.pickle: adjacency list of cluster i after pre-processing.
11. regina_trunc_c[i]_mnodes.pickle: list of manholes in cluster i after pre-processing.
12. regina_trunc_c[i]_mpots.pickle: max potential of each manhole after pre-processing.
13. regina_trunc_c[i]_weights.pickle: min potential of each manhole after pre-processing.

**Python Files description**
1. read_files.py: standalone file and not needed to run the main code. Only used for data-exploration.
2. flow_covering_model_code.py: contains all the functions needed to solve the flow covering model for biosurveillance in wastewater networks. The code is designed to work with two instances: "LA" (Los Angeles) and "regina" (Regina). For "LA", the number of clusters is fixed to 7, and for "regina", the number of clusters is fixed to 10.

**How to Run**

To run the code, you need to modify the orig_task list and provide inst_name. These parameters take the following values:
1. inst_name: either "LA" or "regina".
2. orig_task[0]: Algorithm type (e.g., "M1", "M2", "M3", "M4") - refer to the paper (a brief description also available in the code)
3. orig_task[1]: "flow" or "flow_clus". "flow" solves for the whole network. "flow_clus" solves for a small cluster of the network.
4. orig_task[2]: integer value representing number of sensors to locate
5. orig_task[3]: (Optional) A list of zone indices. This is not optional if "flow_clus" is selected. Zones are 0, 1, 2 for LA, and 1, 2, 3, 4, 5 for Regina.

Example: orig_task = ["M1", "flow", 20], ["M4","flow_clus",100,[2]]

**To save solution map for figures**

Must have done a successful run first and files must be present in results folder. To save solution maps, change the value of save_solution_map to True and run the code with the inst_name and orig_task. At the end, the code will generate two sets of shapefiles in the results folder: one for the solution network, another for the full network it worked on. These shapefiles can be uploaded to ArcGIS software to plot the maps.
