# Standard library imports
from matplotlib.lines import Line2D
import time

# Local imports
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.plotting_module import *
from pomsimulator.modules.text_module import *
from pomsimulator.modules.stats_module import *


def phase_diagram_gen(npz_paths,v_ctt):
    """
        Gathers all concentration arrays, and selects the most predominant species at each pH value.
    Args:
        phase_array_dict: Dictionary of arrays, containing as many concentration array for each species at each pH value
        as initial concentrations calculated and a pH array with all pH values.
        pH: List of pH values
        v_ctt: List of stoichiometries for each molecule

    Returns:
        Bidimensional array containing the indices of the species with most % at each pair of C,pH.
    """
    print("entering phase diagram generation")
    v_ctt_arr_X = np.array([item[0] for item in v_ctt]).reshape(-1, 1)
    v_ctt_arr_M = np.array([item[1] for item in v_ctt]).reshape(-1, 1)

    Ratio_list = list()
    for ii,path in enumerate(npz_paths):

        phase_array_dict = np.load(path)
        labels = phase_array_dict["labels"]
        C_X = phase_array_dict["C_X"]
        C_M = phase_array_dict["C_M"]
        Ratio_list.append(C_M/C_X)
        pH = phase_array_dict["pH"]
        speciation_array = phase_array_dict["SupArray"]
        if ii == 0:
            phase_diagram_X = np.zeros((len(npz_paths), len(pH)), dtype=int)
            phase_diagram_M = np.zeros((len(npz_paths), len(pH)), dtype=int)

        mask_X = get_mask(speciation_array,labels,1.1,C_X,m_idx=0)
        mask_M = get_mask(speciation_array,labels,1.1,C_M,m_idx=1)
        speciation_array_X = speciation_array[:,:,mask_X]
        speciation_array_M = speciation_array[:,:,mask_M]
        Means_X = np.mean(speciation_array_X, axis=2)*v_ctt_arr_X
        Means_M = np.mean(speciation_array_M, axis=2)*v_ctt_arr_M

        phase_diagram_X[ii, :] = np.argmax(Means_X, axis=0)
        phase_diagram_M[ii, :] = np.argmax(Means_M, axis=0)
    return phase_diagram_X,phase_diagram_M,Ratio_list,pH


def main():
    ###########################SETTINGS############################################
    SMALL_SIZE = 10
    MEDIUM_SIZE = 14.5
    BIGGER_SIZE = 14
    plt.rc('font', size=MEDIUM_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    # plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize
    # plt.rc('text', fontsize=MEDIUM_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    Print_logo()
    ###########################PARAMETERS############################################

    formation_labels = Labels_PMo
    speciation_labels = Labels_PMo
    col_dict = None
    v_ctt = [Lab_to_stoich(lab) for lab in speciation_labels]
    npz_info_file = "PMo_npz_info.dat"
    # Load the array with all speciations
    with open(npz_info_file,"r") as infile:
        npz_paths = [line.strip() for line in infile.readlines()]
    phase_diagram_X,phase_diagram_M,Ratio_list,pH = phase_diagram_gen(npz_paths,v_ctt)

    ####################  PLOTTING ######################

    legend = speciation_labels

    if not col_dict:
        nc = len(speciation_labels)
        sample = np.linspace(0,1,nc)
        colors = plt.cm.YlGnBu_r(sample)
        idx_color_dict = {idx:colors[idx] for idx, label in enumerate(speciation_labels)}
    else:
        idx_color_dict = {idx:col_dict[label] for idx, label in enumerate(speciation_labels)}

    rgba_color_dict = {idx: matplotlib.colors.to_rgba(col) for idx, col in idx_color_dict.items()}

    fig,ax = plt.subplots(ncols=2,nrows=1,constrained_layout=True,figsize=(8,4))
    for pp,phase_diagram in enumerate([phase_diagram_X,phase_diagram_M]):
        array_tuple = np.vectorize(rgba_color_dict.get)(phase_diagram)
        tuple_array = np.zeros((len(Ratio_list), len(pH), 4))
        for jj, layer in enumerate(array_tuple):
            tuple_array[:, :, jj] = layer

        color_list = np.unique(phase_diagram).astype(int)
        if col_dict:
            legend_elements = [Line2D([0],[0],marker="o",color=col_dict[speciation_labels[idx]],
                                      label=Lab_to_Formula(legend[idx]),markerfacecolor=col_dict[speciation_labels[idx]],markersize=10)
                               for idx in color_list]
        else:
            legend_elements = [Line2D([0],[0],marker="o",color=colors[idx],
                                      label=Lab_to_Formula(legend[idx]),markerfacecolor=colors[idx],markersize=10)
                               for idx in color_list]


        obj = ax[pp].imshow(tuple_array,
              extent=[max(pH), min(pH), min(Ratio_list), max(Ratio_list)],
              origin='lower',
              aspect='auto', alpha=1,interpolation=None)


        ax[pp].set_xlabel('pH')
        ax[pp].set_ylabel('M/X Ratio')
        ax[pp].set_xlim(min(pH), max(pH))
        leg = ax[pp].legend(handles=legend_elements, loc="center")
        leg.set_draggable(True)

    plt.savefig('../outputs/phase_diagram_PMo.png',dpi=300)

    plt.show()

if __name__ == '__main__':
    main()
