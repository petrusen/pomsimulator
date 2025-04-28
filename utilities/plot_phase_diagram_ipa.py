# Standard library imports
# Local imports
import numpy as np
import matplotlib.pyplot as plt
from pomsimulator.modules.helper_module import phase_diagram_IPA
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.plotting_module import get_color_phase_diagram
from pomsimulator.modules.text_module import Lab_to_stoich,Lab_to_Formula,Print_logo
import re
from configparser import ConfigParser

def main():
    Print_logo()
    config_file = "../inputs/config_W.pomsim"
    config = ConfigParser()
    config.read(config_file)
    ###########################PARAMETERS############################################
    output_path = config["Preparation"]["output_path"]
    system = config["Preparation"]["POM_system"]
    npz_info_file = output_path + "/npz_info_%s.dat" % system
    output_img = "phase_diagram_%s.png" % system

    speciation_labels = config["Speciation"]["speciation_labels"].split(",")
    labels_file = output_path + "/labels_%s.txt" % system
    if speciation_labels[0] == "all":
        with open(labels_file,"r") as flab:
            speciation_labels = [item.strip() for item in flab.readlines()]

    col_dict = None

    v_ctt = [Lab_to_stoich(lab) for lab in speciation_labels]
    # Load the array with all speciations
    with open(npz_info_file,"r") as infile:
        npz_paths = [line.strip() for line in infile.readlines()]
    phase_diagram,C_list,pH = phase_diagram_IPA(npz_paths,v_ctt)

    logC_list = np.log10(C_list)
    ####################  PLOTTING ######################
    color_array,legend_elements = get_color_phase_diagram(phase_diagram,speciation_labels,col_dict)


    # Build the plot
    fig = plt.figure(constrained_layout=True,figsize=(6,4))
    ax = fig.add_subplot()
    obj = ax.imshow(color_array,
          extent=[max(pH), min(pH), min(logC_list), max(logC_list)],
          origin='lower',
          aspect='auto', alpha=1,interpolation=None)


    ax.set_xlabel('pH')
    ax.set_ylabel('$\log_{10}$'+' Concentration (Molar)')
    ax.set_xlim(min(pH), max(pH))
    leg = ax.legend(handles=legend_elements, loc="center")
    leg.set_draggable(True)

    # path for saving the img
    output_img_path = re.sub("array.*npz","",npz_paths[0]) + output_img
    plt.savefig(output_img_path,dpi=300)

    plt.show()

if __name__ == '__main__':
    main()
