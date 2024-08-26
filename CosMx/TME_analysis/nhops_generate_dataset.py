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

## Changed NEIGHBOR_EDGE_CUTOFF to 150 in source
for region_id in region_ids:
    print("Processing %s" % region_id)
    cell_coords_file = os.path.join(raw_data_root, "%s.cell_data.csv" % region_id)
    cell_types_file = os.path.join(raw_data_root, "%s.cell_types.csv" % region_id)
    #cell_biomarker_expression_file = os.path.join(raw_data_root, "%s.expression.csv" % region_id)
    #cell_features_file = os.path.join(raw_data_root, "%s.cell_features.csv" % region_id)
    #voronoi_file = os.path.join(raw_data_root, "%s.json" % region_id)
    voronoi_polygon_img_output = os.path.join(fig_save_root, "%s_voronoi.png" % region_id)
    graph_img_output = os.path.join(fig_save_root, "%s_graph.png" % region_id)
    graph_output = os.path.join(nx_graph_root, "%s.gpkl" % region_id)
    mapping_output = os.path.join(map_save_root, "%s.csv" % region_id)
    if not os.path.exists(mapping_output):
        G,mapping = spacegm.construct_graph_for_region(
            region_id,
            cell_coords_file=cell_coords_file,
            cell_types_file=cell_types_file,
            #cell_biomarker_expression_file=cell_biomarker_expression_file,
            #cell_features_file=cell_features_file,
            #voronoi_file=voronoi_file,
            graph_output=graph_output,
            voronoi_polygon_img_output=voronoi_polygon_img_output,
            graph_img_output=graph_img_output,
            figsize=10,
            return_cell_mapping = True)
        with open(mapping_output, 'w') as csv_file:  
            writer = csv.writer(csv_file)
            for key, value in mapping.items():
                writer.writerow([key, value])    
plt.close()
print("Finished!")

