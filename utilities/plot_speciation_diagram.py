# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from pomsimulator.modules.plotting_module import plot_speciation
from pomsimulator.modules.helper_module import load_array,get_C0
from pomsimulator.modules.text_module import Print_logo, Lab_to_Formula
from configparser import ConfigParser

def main():
    Print_logo()
    config_file = "../inputs/config_PMo.pomsim"
    config = ConfigParser()
    config.read(config_file)
    output_path = config["Preparation"]["output_path"]
    system = config["Preparation"]["POM_system"]
    npz_file = output_path + "/" + config["Speciation"]["npz_file"]
    m_idx = int(config["Speciation"]["m_idx"])
    output_img = output_path + "/Speciation_Diagram_%s.png" % system

    fig,ax = plt.subplots(1,1,figsize=(6.5,4),constrained_layout=True)

    conc_arr,index_arr,C_ref,pH,speciation_labels = load_array(npz_file)
    C0 = get_C0(C_ref,m_idx=m_idx)
    conc_means = np.mean(conc_arr,axis=2)
    plot_speciation(conc_means,speciation_labels,pH,C0,None,ax=ax,m_idx=m_idx)

    handle,labels = ax.get_legend_handles_labels()
    ax.legend(handle,[Lab_to_Formula(lab) for lab in labels])

    plt.savefig(output_img,dpi=300)
    plt.show()

if __name__ == '__main__':
    main()
