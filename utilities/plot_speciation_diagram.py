# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from pomsimulator.modules.DataBase import Labels_W_good
from pomsimulator.modules.plotting_module import plot_speciation
from pomsimulator.modules.stats_module import load_array
from pomsimulator.modules.text_module import Print_logo, Lab_to_Formula

def main():
    Print_logo()
    path = "../outputs/W_Array.npz"
    speciation_labels = Labels_W_good

    fig,ax = plt.subplots(1,1,figsize=(6.5,4),constrained_layout=True)

    conc_arr,index_arr,pH = load_array(path)

    conc_means = np.mean(conc_arr,axis=2)

    plot_speciation(conc_means,speciation_labels,pH,0.1,None,ax=ax)

    handle,labels = ax.get_legend_handles_labels()
    ax.legend(handle,[Lab_to_Formula(lab) for lab in labels])

    plt.savefig('../outputs/Speciation_diagram_W.png',dpi=300)
    plt.show()

if __name__ == '__main__':
    main()
