#Standard library imports
import sys
import os
import re
import time
import pkg_resources as pkgr
from configparser import ConfigParser
#Third-party imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#Local imports
from pomsimulator.modules.stats_module import *
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.text_module import Print_logo, Lab_to_Formula
from pomsimulator.modules.helper_module import load_array,get_C0,get_config_to_dict
from pomsimulator.modules.plotting_module import plot_cluster_means,plot_speciation


def SM_clustering(config_dict):
    """
    Performs clustering on the data.

    Args:
         config_dict: Configuration dictionary containing parameters.
    Returns:
        str: Status message
    """

    try:

        Print_logo()
        run_by_gui = True

        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        ############## INPUT VARIABLES ######################

        output_path = pkgr.resource_filename(__name__,config_dict["Preparation"]["output_path"])
        system = config_dict["Preparation"]["POM_system"]

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        path_to_npz = output_path + "/" + config_dict["Clustering"]["npz_cluster_file"]
        path_to_features = output_path + "/" + config_dict["Clustering"]["features_file"]

        clust_dir = output_path + "/" + config_dict["Clustering"]["cluster_dir"]
        clust_path = clust_dir

        n_clusters = int(float(config_dict["Clustering"]["n_clusters"]))
        normalize = config_dict["Clustering"]["normalize_feats"]
        feats_list = config_dict["Clustering"]["feats_list"].split(",")
        col_dict = None
        if config_dict["Visualization"]["col_dict"] != "":
            col_dict_name = config_dict["Visualization"]["col_dict"]
            col_dict = color_dictionaries[col_dict_name]


        # 1) Load Array ############################################################################################
        print("1) Loading Arrays")

        t0 = time.time()
        SuperArr,IndexArr,C_ref,pH,labels = load_array(path_to_npz)

        if config_dict["Visualization"]["plot_list"]!= "all":
            plot_list = config_dict["Visualization"]["plot_list"].split(",")
        else:
            plot_list = list(labels)

        t1 = time.time()
        print("Array loaded in %.1f s" % (t1 - t0))

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        mask = np.full(SuperArr.shape[2],True)
        for ii in range(len(system.split("_"))):
            C0 = get_C0(C_ref,ii)
            mask *= mask_models(SuperArr,labels,threshold=1.1,C=C0,m_idx=ii)
        SuperArr = SuperArr[:,:,mask]
        IndexArr = IndexArr[IndexArr != -1][mask]

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 2) Computing features
        print("2) Computing features")

        t0 = time.time()
        if os.path.exists(path_to_features):
            df_full = pd.read_csv(path_to_features, index_col=0)
            print("Read existing features (%s), %d entries" % (path_to_features,df_full.shape[0]))
        else:
            df_full = get_features_array(SuperArr, labels, pH)
            df_full.loc[:,"mod_idx"] = IndexArr
            df_full.to_csv(path_to_features)

        col_sel = [col for col in df_full.columns if col.split("-")[0] in feats_list]
        areas = [col for col in df_full.columns if "area" in col]
        df = df_full.loc[:, col_sel]

        t1 = time.time()
        print("Featurized after %.1f s" % (t1-t0))

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 3) Clustering
        print("3) Clustering")

        if not os.path.exists(clust_path):
            os.makedirs(clust_path)
        else:
            counter = 1
            while os.path.exists(clust_path):
                clust_path = clust_dir + ".%03d" % counter
                counter += 1
            else:
                os.makedirs(clust_path)

        time0 = time.time()
        cluster_info = get_clusters(df,n_clusters,normalize=normalize)

        print(["%.2f" % val for val in cluster_info["PCA"].explained_variance_ratio_])
        print("Generating groups")
        groups = get_cluster_members(cluster_info)
        print("Clusterised after %.2f"%(time.time()-time0))

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        with open(clust_path + "/groups.csv", 'w') as fout:
            for ii, grp in enumerate(groups):
                fout.write(str(ii)+','+','.join([str(gr) for gr in grp]) + "\n")
                print("Cluster %d, %d elements" % (ii, len(grp)))

        for ii in range(len(system.split("_"))):
            C0 = get_C0(C_ref,ii)
            fig,axd = plot_cluster_means(SuperArr,groups,labels,pH,C0,
                                        col_dict=col_dict,plot_list=plot_list,
                                        target_shape=None,m_idx=ii)

            plt.savefig(clust_path + '/clusters_speciation_idx%d.svg' % ii ,dpi=300,transparent=False)

        return("Normal termination")

    except Exception as e:
        print("An error occurred during the clustering: ", str(e))
        print(f"Error {e}")
        return "Error"

def cluster_selection(config_dict):
    """
    Selects clusters for further analysis.

    Args:
         config_dict: Configuration dictionary containing parameters.
         """
    try:
        Print_logo()
        run_by_gui = True

        if type(config_dict)!= dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        output_path = pkgr.resource_filename(__name__,config_dict["Preparation"]["output_path"])
        system = config_dict["Preparation"]["POM_system"]

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        clust_dir = output_path + "/" + config_dict["Clustering"]["cluster_dir"]
        path_to_npz = output_path + "/" + config_dict["Clustering"]["npz_cluster_file"]

        groups_file_path = clust_dir + "/groups.csv"
        sel_groups_idx = [int(float(item)) for item in config_dict["Clustering"]["sel_groups"].split(",")]
        features_file_path = output_path + "/" + config_dict["Clustering"]["features_file"]

        SuperArr,IndexArr,C_ref,pH,labels = load_array(path_to_npz)
        ### Apply masks
        mask = np.full(SuperArr.shape[2],True)
        for ii in range(len(system.split("_"))):
            C0 = get_C0(C_ref,ii)
            mask *= mask_models(SuperArr,labels,threshold=1.1,C=C0,m_idx=ii)
        SuperArr = SuperArr[:,:,mask]
        IndexArr = IndexArr[IndexArr != -1][mask]

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        Groups_dict = dict()
        with open(groups_file_path,'r') as infile:
            lines = [line.strip().split(',') for line in infile.readlines()]
            for line in lines:
                Groups_dict[int(line[0])]= [int(element) for element in line[1:]]

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        group_indices = sum([Groups_dict[jj] for jj in sel_groups_idx],[])
        NewArr = SuperArr[:,:,np.array(group_indices)]
        NewIndArr = IndexArr[np.array(group_indices)]

        np.savez_compressed(clust_dir + '/NewArr.npz', SupArray=NewArr, IndexArray=NewIndArr,
                            pH=pH, labels=labels, selected_groups=sel_groups_idx,
                            original_npz=path_to_npz,C_ref=C_ref)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        if features_file_path:
            df_full = pd.read_csv(features_file_path, index_col=0)
            df_new = df_full.iloc[group_indices,:].copy()
            df_new.to_csv(clust_dir + "/new_features.csv")

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        model_idxs = [IndexArr[k] for k in sorted(group_indices)]
        with open(clust_dir + "/sel_model_indices.pomsim","w") as fmod:
            fmod.write("\n".join([str(idx) for idx in model_idxs]))

        print("Selected model indices saved to ", clust_dir)

        return("Normal termination")

    except Exception as e:
        print("An error occurred during the cluster selection: ", str(e))
        print(f"Error {e}")
        return "Error"

def clust_filtering(config_dict):
    """
    Performs filtering on the selected clusters.

    Args:
         config_dict: Configuration dictionary containing parameters.
         """
    try:
        Print_logo()
        run_by_gui = True
        if type(config_dict)!= dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        output_path = pkgr.resource_filename(__name__,config_dict["Preparation"]["output_path"])
        system = config_dict["Preparation"]["POM_system"]

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        clust_dir = output_path + "/" + config_dict["Clustering"]["cluster_dir"]
        path_to_npz = clust_dir + "/NewArr.npz"

        SuperArr,IndexArr,C_ref,pH,labels = load_array(path_to_npz)
        labels = list(labels)
        filter_path = clust_dir + "/filtering"
        m_idx = int(config_dict["Speciation"]["m_idx"])

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        col_dict = None
        if config_dict["Visualization"]["col_dict"] != "":
            col_dict_name = config_dict["Visualization"]["col_dict"]
            col_dict = color_dictionaries[col_dict_name]

        if config_dict["Visualization"]["plot_list"]!= "all":
            plot_list = config_dict["Visualization"]["plot_list"].split(",")
        else:
            plot_list = list(labels)

        if config_dict["Visualization"]["boxplot_list"]!= "all":
            boxplot_list = config_dict["Visualization"]["boxplot_list"].split(",")
        else:
            boxplot_list = labels

        C0 = get_C0(C_ref,m_idx)
        os.makedirs(filter_path,exist_ok=True)

        for sel_spc in boxplot_list:
            if sel_spc not in labels:
                print(sel_spc, "was not found")
                continue

            if stop_signal and os.path.exists(stop_signal):
                print("Function stopped by user")
                return "Stopped"

            sel_ndx = labels.index(sel_spc)
            FilterSelArr = boxplot_filtering(sel_ndx, pH, labels, SuperArr)
            print(sel_spc,FilterSelArr.shape)
            fig = plt.figure(constrained_layout=True,figsize=(7,7))

            ax=fig.subplot_mosaic([[0],[1]],gridspec_kw={"height_ratios":[1,0.35]})

            plot_speciation(np.mean(FilterSelArr,axis=2),labels,pH,C0,
                            plot_list=plot_list,ax=ax[0],m_idx=m_idx,err_arr=np.std(FilterSelArr,axis=2),col_dict=col_dict)

            ax[0].set_title("%s"%sel_spc)

            handles, leg_labels = ax[0].get_legend_handles_labels()
            leg_labels = [Lab_to_Formula(lab) for lab in leg_labels]
            ax[1].legend(handles,leg_labels,ncols=4,loc='center',fontsize=9)
            ax[1].axis("off")
            # remove the underlying axes for lower part

            plt.savefig(filter_path + "/filt_by_%s.svg" % sel_spc,
                        dpi=300)

        return "Normal termination"

    except Exception as e:
        print("An error occurred during the filtering: ", str(e))
        print(f"Error {e}")
        return "Error"