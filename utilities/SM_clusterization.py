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
from pomsimulator.modules.text_module import Print_logo
from pomsimulator.modules.helper_module import load_array,get_C0
from pomsimulator.modules.plotting_module import plot_cluster_means


def main():
    Print_logo()

    ############## INPUT VARIABLES ######################
    config_file = pkgr.resource_filename(__name__, "../inputs/config_PMo.pomsim")
    config = ConfigParser()
    config.read(config_file)

    output_path = pkgr.resource_filename(__name__, config["Preparation"]["output_path"])
    system = config["Preparation"]["POM_system"]

    path_to_npz = output_path + "/" + config["Clustering"]["npz_cluster_file"]
    path_to_features = output_path + "/" + config["Clustering"]["features_file"]

    clust_dir = output_path + "/" + config["Clustering"]["cluster_dir"]
    clust_path = clust_dir

    n_clusters = int(config["Clustering"]["n_clusters"])
    normalize = config["Clustering"].getboolean("normalize_feats")
    feats_list = config["Clustering"]["feats_list"].split(",")

    col_dict = None
    plot_list = None

    # 1) Load Array
    print("1) Loading Arrays")

    t0 = time.time()
    SuperArr,IndexArr,C_ref,pH,labels = load_array(path_to_npz)

    t1 = time.time()
    print("Array loaded in %.1f s" % (t1 - t0))

    mask = np.full(SuperArr.shape[2],True)
    for ii in range(len(system.split("_"))):
        C0 = get_C0(C_ref,ii)
        mask *= mask_models(SuperArr,labels,threshold=1.1,C=C0,m_idx=ii)
    SuperArr = SuperArr[:,:,mask]
    IndexArr = IndexArr[IndexArr != -1][mask]

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

    with open(clust_path + "/groups.csv", 'w') as fout:
        for ii, grp in enumerate(groups):
            fout.write(str(ii)+','+','.join([str(gr) for gr in grp]) + "\n")
            print("Cluster %d, %d elements" % (ii, len(grp)))

    if not plot_list:
        plot_list = list(labels)

    for ii in range(len(system.split("_"))):
        C0 = get_C0(C_ref,ii)
        fig,axd = plot_cluster_means(SuperArr,groups,labels,pH,C0,
                                    col_dict=col_dict,plot_list=plot_list,
                                    target_shape=None,m_idx=ii)
        plt.savefig(clust_path + '/clusters_speciation_idx%d.svg' % ii ,dpi=300,transparent=False)
    print("Normal termination")
if __name__ == "__main__":
    main()
