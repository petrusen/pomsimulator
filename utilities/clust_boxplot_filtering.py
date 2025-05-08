#Standard library imports
import os
from configparser import ConfigParser
import pkg_resources as pkgr
#Third-party imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Local imports
from pomsimulator.modules.DataBase import Col_Dict_PMo
from pomsimulator.modules.text_module import Print_logo,Lab_to_Formula
from pomsimulator.modules.stats_module import *
from pomsimulator.modules.helper_module import load_array,get_C0
from pomsimulator.modules.plotting_module import plot_speciation



def main():
    Print_logo()
    config_file = pkgr.resource_filename(__name__, "../inputs/config_PMo.pomsim")
    config = ConfigParser()
    config.read(config_file)

    output_path = pkgr.resource_filename(__name__, config["Preparation"]["output_path"])
    system = config["Preparation"]["POM_system"]
    clust_dir = output_path + "/" + config["Clustering"]["cluster_dir"]
    path_to_npz = clust_dir + "/NewArr.npz"

    SuperArr,IndexArr,C_ref,pH,labels = load_array(path_to_npz)
    labels = list(labels)
    filter_path = clust_dir + "/filtering"
    m_idx = int(config["Speciation"]["m_idx"])

    col_dict = None
    plot_list = None
    boxplot_list = None
    C0 = get_C0(C_ref,m_idx)
    os.makedirs(filter_path,exist_ok=True)

    if not boxplot_list:
        boxplot_list = labels

    if not plot_list:
        plot_list = labels

    for sel_spc in boxplot_list:
        if sel_spc not in labels:
            print(sel_spc, "was not found")
            continue

        sel_ndx = labels.index(sel_spc)
        FilterSelArr = boxplot_filtering(sel_ndx, pH, labels, SuperArr)
        print(sel_spc,FilterSelArr.shape)
        fig = plt.figure(constrained_layout=True,figsize=(7,7))

        ax=fig.subplot_mosaic([[0],[1]],gridspec_kw={"height_ratios":[1,0.35]})

        plot_speciation(np.mean(FilterSelArr,axis=2),labels,pH,C0,
                        plot_list=None,ax=ax[0],m_idx=m_idx,err_arr=np.std(FilterSelArr,axis=2),col_dict=col_dict)
        ax[0].set_title("%s"%sel_spc)

        handles, leg_labels = ax[0].get_legend_handles_labels()
        leg_labels = [Lab_to_Formula(lab) for lab in leg_labels]
        ax[1].legend(handles,leg_labels,ncols=4,loc='center',fontsize=9)
        ax[1].axis("off")
        # remove the underlying axes for lower part

        plt.savefig(filter_path + "/filt_by_%s.svg" % sel_spc,
                    dpi=300)
    print("Normal termination")
if __name__ == "__main__":
    main()
