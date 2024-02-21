# Standard library imports
from matplotlib.lines import Line2D
import time

# Local imports
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.plotting_module import *
from pomsimulator.modules.text_module import *
from pomsimulator.modules.stats_module import *


def phase_diagram_gen(phase_array_dict,pH,v_ctt):
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
    C_list = phase_array_dict["conc_array"]
    phase_diagram = np.zeros((len(C_list), len(pH)), dtype=int)
    v_ctt_arr = np.array([item[0] for item in v_ctt]).reshape(-1,1)
    for i, C in enumerate(C_list):
        init_time = time.time()
        speciation_array = phase_array_dict["array_%d" % i]
        half = time.time()
        mask = get_mask(speciation_array,1.1*C)
        end = time.time()
        speciation_array = speciation_array[:,:,mask]
        Means = np.mean(speciation_array, axis=2)
        Means = Means*v_ctt_arr
        phase_ratio = np.argmax(Means, axis=0)
        phase_diagram[i, :] = phase_ratio
    return phase_diagram


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

    formation_labels = Labels_W
    speciation_labels = Labels_W_good
    col_dict = None
    v_ctt = [Lab_to_stoich(lab) for lab in speciation_labels]
    # Load the array with all speciations
    phase_array = np.load("../outputs/phase_array.npz")
    pH = phase_array["pH"]
    C_list = phase_array["conc_array"]
    phase_diagram = phase_diagram_gen(phase_array,pH,v_ctt)

    logC_list = np.log10(C_list)
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

    array_tuple = np.vectorize(rgba_color_dict.get)(phase_diagram)
    tuple_array = np.zeros((len(C_list), len(pH), 4))
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

    fig = plt.figure(constrained_layout=True,figsize=(6,4))
    ax = fig.add_subplot()
    obj = ax.imshow(tuple_array,
          extent=[max(pH), min(pH), min(logC_list), max(logC_list)],
          origin='lower',
          aspect='auto', alpha=1,interpolation=None)


    ax.set_xlabel('pH')
    ax.set_ylabel('$\log_{10}$'+' Concentration (Molar)')
    ax.set_xlim(min(pH), max(pH))
    leg = ax.legend(handles=legend_elements, loc="center")
    leg.set_draggable(True)

    plt.savefig('../outputs/phase_diagram_W.png',dpi=300)

    plt.show()

if __name__ == '__main__':
    main()
