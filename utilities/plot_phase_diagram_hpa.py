#Standard library imports
import re
from configparser import ConfigParser
import pkg_resources as pkgr
#Third-party imports
import matplotlib.pyplot as plt
# Local imports
from pomsimulator.modules.helper_module import phase_diagram_HPA
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.plotting_module import get_color_phase_diagram
from pomsimulator.modules.text_module import Lab_to_stoich,Print_logo


def main():
    Print_logo()
    config_file = pkgr.resource_filename(__name__, "../inputs/config_PMo.pomsim")
    config = ConfigParser()
    config.read(config_file)
    ###########################PARAMETERS############################################

    output_path = pkgr.resource_filename(__name__, config["Preparation"]["output_path"])
    system = config["Preparation"]["POM_system"]
    npz_info_file = output_path + "/npz_info_%s.dat" % system

    output_img = "phase_diagram_%s.png" % system

    speciation_labels = config["Speciation"]["speciation_labels"].split(",")
    labels_file = output_path + "/labels_%s.txt" % system
    if speciation_labels[0] == "all":
        with open(labels_file,"r") as flab:
            speciation_labels = [item.strip() for item in flab.readlines()]

    col_dict = Col_Dict_PMo

    v_ctt = [Lab_to_stoich(lab) for lab in speciation_labels]
    # Load the array with all speciations
    with open(npz_info_file,"r") as infile:
        npz_paths = [line.strip() for line in infile.readlines()]

    phase_diagram_X,phase_diagram_M,Ratio_list,pH = phase_diagram_HPA(npz_paths,v_ctt)

    ####################  PLOTTING ######################


    fig,ax = plt.subplots(ncols=2,nrows=1,constrained_layout=True,figsize=(8,4))
    for pp,phase_diagram in enumerate([phase_diagram_X,phase_diagram_M]):
        color_array,legend_elements = get_color_phase_diagram(phase_diagram,speciation_labels,col_dict)
        obj = ax[pp].imshow(color_array,
              extent=[max(pH), min(pH), min(Ratio_list), max(Ratio_list)],
              origin='lower',
              aspect='auto', alpha=1,interpolation=None)


        ax[pp].set_xlabel('pH')
        ax[pp].set_ylabel('M/X Ratio')
        ax[pp].set_xlim(min(pH), max(pH))
        leg = ax[pp].legend(handles=legend_elements, loc="center")
        leg.set_draggable(True)

    # path for saving the img

    output_img_path = re.sub("array.*npz","",npz_paths[0]) + ["X_","M_"][pp] + output_img
    plt.savefig(output_img_path,dpi=300)

    plt.show()
    print("Normal termination")
if __name__ == '__main__':
    main()
