import numpy as np
import pandas as pd

from pomsimulator.modules.text_module import Print_logo
from pomsimulator.modules.stats_module import *
from pomsimulator.modules.helper_module import load_array

output_path = "../outputs/PMo_data/PMo_clusterization"
path_to_npz = output_path + "/../PMo_Array.npz"
groups_file_path = output_path + "/groups.csv"
sel_groups_idx = [0,1,2]
features_file_path = output_path + "/../features_PMo.csv"

SuperArr,IndexArr,C_ref,pH,labels = load_array(path_to_npz)

Groups_dict = dict()
with open(groups_file_path,'r') as infile:
    lines = [line.strip().split(',') for line in infile.readlines()]
    for line in lines:
        Groups_dict[int(line[0])]= [int(element) for element in line[1:]]


group_indices = sum([Groups_dict[jj] for jj in sel_groups_idx],[])
NewArr = SuperArr[:,:,np.array(group_indices)]
NewIndArr = IndexArr[np.array(group_indices)]

np.savez_compressed(output_path + '/NewArr.npz', SupArray=NewArr, IndexArray=NewIndArr,
                    pH=pH, labels=labels, selected_groups=sel_groups_idx,
                    original_npz=path_to_npz,C_ref=C_ref)

if features_file_path:
    df_full = pd.read_csv(features_file_path, index_col=0)
    df_new = df_full.iloc[group_indices,:].copy()
    df_new.to_csv(output_path + "/new_features.csv")