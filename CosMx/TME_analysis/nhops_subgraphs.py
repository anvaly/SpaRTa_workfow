#!/usr/bin/env python
# coding: utf-8

#imports
import spacegm
import argparse

import warnings
warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import torch.nn as nn
import torch
import umap
import scipy
import csv
import seaborn as sns
plt.rcParams['figure.figsize'] = [15, 9]


parser = argparse.ArgumentParser()
parser.add_argument("-n", "--nhops", help="number of hops", type=int)
parser.add_argument("-s", "--nsize", help="neighborhood size", type=int)
args = parser.parse_args()

nhops=args.nhops
nsize=args.nsize

##SPACE-GM dataset object
raw_data_root = "/opt/dlami/nvme/cosmx/data/intermediate/TME_analysis/tme_subtype/"
dataset_root = "/opt/dlami/nvme/cosmx/data/intermediate/TME_analysis/tme_subtype/dataset/"
nx_graph_root = os.path.join(dataset_root, "graph")
fig_save_root = os.path.join(dataset_root, "fig")
map_save_root = os.path.join(dataset_root, "maps")
model_save_root = os.path.join(dataset_root, 'model')

os.makedirs(nx_graph_root, exist_ok=True)
os.makedirs(fig_save_root, exist_ok=True)
os.makedirs(model_save_root, exist_ok=True)
os.makedirs(map_save_root, exist_ok=True)

region_ids = [
    'slide1_FOV1',
    'slide1_FOV2',
    'slide1_FOV3',
    'slide1_FOV4',
    'slide1_FOV5',
    'slide1_FOV6',
    'slide1_FOV7',
    'slide1_FOV8',
    'slide1_FOV9',
    'slide1_FOV10',
    'slide1_FOV11',
    'slide1_FOV12',
    'slide1_FOV13',
    'slide1_FOV14',
    'slide1_FOV15',
    'slide1_FOV16',
    'slide1_FOV17',
    'slide1_FOV18',
    'slide1_FOV19',
    'slide1_FOV20',
    'slide1_FOV21',
    'slide1_FOV22',
    'slide1_FOV23',
    'slide2_FOV1',
    'slide2_FOV2',
    'slide2_FOV3',
    'slide2_FOV4',
    'slide2_FOV5',
    'slide2_FOV6',
    'slide2_FOV7',
    'slide2_FOV8',
    'slide2_FOV9',
    'slide2_FOV10',
    'slide2_FOV11',
    'slide2_FOV12',
    'slide2_FOV13',
    'slide2_FOV14',
    'slide2_FOV15',
    'slide2_FOV16',
    'slide2_FOV17',
    'slide2_FOV18',
    'slide2_FOV19',
    'slide2_FOV20',
    'slide2_FOV21',
    'slide3_FOV1',
    'slide3_FOV2',
    'slide3_FOV3',
    'slide3_FOV4',
    'slide3_FOV5',
    'slide3_FOV6',
    'slide3_FOV7',
    'slide3_FOV8',
    'slide3_FOV9',
    'slide3_FOV10',
    'slide3_FOV11',
    'slide3_FOV12',
    'slide3_FOV13',
    'slide3_FOV14',
    'slide3_FOV15',
    'slide3_FOV16',
    'slide3_FOV17',
    'slide3_FOV18',
    'slide3_FOV19',
    'slide3_FOV20',
    'slide4_FOV1',
    'slide4_FOV2',
    'slide4_FOV3',
    'slide4_FOV4',
    'slide4_FOV5',
    'slide4_FOV6',
    'slide4_FOV7',
    'slide4_FOV8',
    'slide4_FOV9',
    'slide4_FOV10',
    'slide4_FOV11',
    'slide4_FOV12',
    'slide4_FOV13',
    'slide4_FOV14',
    'slide4_FOV15',
    'slide4_FOV16',
    'slide4_FOV17',
    'slide4_FOV18',
    'slide4_FOV19',
    'slide4_FOV20',
    'slide4_FOV21',
    'slide5_FOV14',
    'slide5_FOV15',
    'slide5_FOV16',
    'slide5_FOV17',
    'slide5_FOV18',
    'slide5_FOV19',
    'slide5_FOV20',
    'slide5_FOV21',
    'slide5_FOV22',
    'slide5_FOV23',
    'slide5_FOV24',
    'slide5_FOV25',
    'slide5_FOV26',
    'slide5b2_FOV11',
    'slide5b2_FOV12',
    'slide5b2_FOV13',
    'slide5b2_FOV14',
    'slide5b2_FOV15',
    'slide5b2_FOV16',
    'slide5b2_FOV17',
    'slide5b2_FOV18',
    'slide5b2_FOV19',
    'slide5b2_FOV20',
    'slide5b2_FOV21',
    'slide5b2_FOV22',
    'slide5b2_FOV23',
    'slide5b2_FOV24',
    'slide5b3_FOV1',
    'slide5b3_FOV2',
    'slide5b3_FOV3',
    'slide5b3_FOV6',
    'slide5b3_FOV7',
    'slide5b3_FOV8',
    'slide5b3_FOV10',
    'slide5b3_FOV11',
    'slide5b3_FOV13',
    'slide5b3_FOV14',
    'slide5b3_FOV15',
    'slide5b3_FOV16',
    'slide5b3_FOV18',
    'slide5b3_FOV19',
    'slide5b3_FOV20',
    'slide6_FOV1',
    'slide6_FOV2',
    'slide6_FOV3',
    'slide6_FOV4',
    'slide6_FOV5',
    'slide6_FOV6',
    'slide6_FOV7',
    'slide6_FOV8',
    'slide6_FOV9',
    'slide6_FOV10',
    'slide6_FOV11',
    'slide6_FOV12',
    'slide6_FOV13',
    'slide6_FOV14',
    'slide6_FOV15',
    'slide6_FOV16',
    'slide6_FOV17',
    'slide6_FOV18',
    'slide6_FOV19',
    'slide6_FOV20',
]

dataset_kwargs = {
    'transform': [],
    'pre_transform': None,
    'raw_folder_name': 'graph',  # os.path.join(dataset_root, "graph") is the folder where we saved nx graphs
    'processed_folder_name': 'tg_graph',  # processed dataset files will be stored here
    'node_features': ["cell_type", "neighborhood_composition", "center_coord"],  # There are all the cellular features that we want the dataset to compute
    'edge_features': ["edge_type", "distance"],  # edge (cell pair) features
    'subgraph_size': nhops,  # how many hops? 
    'subgraph_source': 'on-the-fly',
    'subgraph_allow_distant_edge': True,
    'subgraph_radius_limit': 200.,
}

feature_kwargs = {
#    "biomarker_expression_process_method": "linear",
#    "biomarker_expression_lower_bound": 0,
#    "biomarker_expression_upper_bound": 18,
    "neighborhood_size": nsize
}
dataset_kwargs.update(feature_kwargs)


dataset = spacegm.CellularGraphDataset(dataset_root, **dataset_kwargs)
dataset.cell_type_mapping

def get_subgraph_ids(sg, fov):
    mapping_file = os.path.join(map_save_root, "%s.csv" % sg.region_id)
    with open(mapping_file) as csv_file:
        reader = csv.reader(csv_file)
        mapping = dict(reader)
    sg_nodes = dataset.calculate_subgraph_nodes(fov_map[sg.region_id], sg.original_center_node)
    return [mapping[str(val)] for val in sg_nodes.tolist()]


#get number of cells/TMEs per FOV
num_subgraphs_by_fov = {}
fov_map = {}
for fov in range(len(region_ids)):
    num_subgraphs_by_fov[fov] = dataset.get_full(fov).x.shape[0]
    fov_map[dataset.get_full(fov).region_id] = fov
    
#threshold for excluding TMEs
exclude_threshold = 0.2

#define data structure
comp_dataset = {
    "x_normalized": {}, 
    "fov": {},
    "center_cell_identity": {},
    "original_sg_idx": {},
    "cell_members": {},
    "center_cell": {}
}

num_cell_types = len(dataset.cell_type_mapping.keys())

#iterate over FOVs and their subgraphs
for fov, num_sg in num_subgraphs_by_fov.items():
    comp_dataset["x_normalized"][fov] = []
    comp_dataset["fov"][fov] = []
    comp_dataset["center_cell_identity"][fov] = []
    comp_dataset["original_sg_idx"][fov] = []
    comp_dataset["cell_members"][fov] = []
    comp_dataset["center_cell"][fov] = []
    for sg_idx in range(num_sg):
        sg = dataset.get_subgraph(fov, sg_idx)
        #count elements
        cell_list = [sg.x[cell,0].item()  for cell in range(sg.x.shape[0]) if (cell != sg["center_node_index"])  ]
        cell_count_dict = Counter(cell_list)
        #create composition vector
        composition_vector = [cell_count_dict[cell] / (sg.x.shape[0]-1) if cell in cell_count_dict.keys() else 0 for cell in range(num_cell_types)]    
        #quality control
        if (composition_vector[8] + composition_vector[11] < exclude_threshold):            
            #append to dataset obj
            comp_dataset["x_normalized"][fov].append(composition_vector)
            comp_dataset["fov"][fov].append(sg.region_id)
            comp_dataset["center_cell_identity"][fov].append(sg.x[sg.center_node_index, 0].item())
            comp_dataset["original_sg_idx"][fov].append(sg_idx)
            comp_dataset["cell_members"][fov].append(get_subgraph_ids(sg, fov))
            comp_dataset["center_cell"][fov].append(get_subgraph_ids(sg, fov)[sg.center_node_index])
    #x into sparse array
    #comp_dataset["x_normalized"][fov] = np.array(comp_dataset["x_normalized"][fov])
    #comp_dataset["x_normalized"][fov] =  scipy.sparse.csr_array(comp_dataset["x_normalized"][fov])
    print("Processed FOV:", fov)  


##Building Composition Matrix
total_subgraphs = 0
num_valid_tmes_fov = {}

#describe dataset
for fov in range(len(region_ids)):    
    num_valid_tmes_fov[fov] = len(comp_dataset["original_sg_idx"][fov])
    total_subgraphs += num_valid_tmes_fov[fov]

total_subgraphs

#build x matrix
for fov in num_valid_tmes_fov.keys():
    if len(comp_dataset["x_normalized"][fov]) != 0:
        if "x" not in locals():
            x = np.array(comp_dataset["x_normalized"][fov])
            #retain original index in dataset object
            x_original_sg = np.array([(fov, orig_idx) for orig_idx in comp_dataset["original_sg_idx"][fov]])
        else:
            x = np.concatenate((x, np.array(comp_dataset["x_normalized"][fov])))
            x_original_sg = np.concatenate((x_original_sg, np.array([(fov, orig_idx) for orig_idx in comp_dataset["original_sg_idx"][fov]])))
 

#imports
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

plt.rcParams['figure.figsize'] = [15, 9]
dataset.cell_type_mapping

#save cluster memberships for each k
k_cluster_memberships = {}

#cluster all patients
ks=[4, 8, 15, 20]
for k in ks:
    kmeans = KMeans(n_clusters=k, random_state=123)
    cluster_mem = kmeans.fit_predict(x, sample_weight=None)

    #save memberships
    k_cluster_memberships[k] = cluster_mem

    #get centroids
    centroids = kmeans.cluster_centers_

    #visualize centroids as heatmap
    cmap = sns.color_palette("rocket_r", as_cmap=True)
    sns.heatmap(centroids, cmap=cmap, annot=True, fmt=".2f")
    plt.xlabel("Cell_types")
    plt.xticks(range(centroids.shape[1]), [cell_type for cell_type in dataset.cell_type_mapping.keys()], rotation=45)
    plt.ylabel("Clusters")
    plt.title(f"Heatmap of cluster centroids with k = {k} based on {x.shape[0]} random TMEs")
    plt.savefig(f"/opt/dlami/nvme/cosmx/data/intermediate/TME_analysis/tme_subtype/outputs/composition_full_heatmap_{k}_{nhops}hop.png")
    plt.show()


#save comp_full_qc_20
x_df = pd.DataFrame(x, columns= [cell_type for cell_type in dataset.cell_type_mapping.keys()])
x_inds_df = pd.DataFrame(x_original_sg, columns=["FOV", "Subgraph ID"])
x_memberships = pd.DataFrame(k_cluster_memberships)
x_memberships.columns = ["k = " + str(k) for k in k_cluster_memberships.keys()]
x_df = pd.concat((x_inds_df, x_memberships, x_df), axis=1)

x_df.to_csv("/opt/dlami/nvme/cosmx/data/intermediate/TME_analysis/tme_subtype/composition_full_qc_20_" + str(nhops) + "hop_" + str(nsize) + "nsize.csv")


##Collect metadata
pds = []
for fov in num_valid_tmes_fov.keys():
    if len(comp_dataset["x_normalized"][fov]) != 0:
        pds.append(pd.DataFrame({"cell_members" : comp_dataset["cell_members"][fov],
                                 "fov" : comp_dataset["fov"][fov],
                                 "center_cell" : comp_dataset["center_cell"][fov],
                                 "FOV" : fov, 
                                 "Subgraph ID" : comp_dataset["original_sg_idx"][fov]}))
        
meta = pd.concat(pds)
meta.to_csv("/opt/dlami/nvme/cosmx/data/intermediate/TME_analysis/tme_subtype/composition_meta_" + str(nhops) + "hop_" + str(nsize) + "nsize.csv")
