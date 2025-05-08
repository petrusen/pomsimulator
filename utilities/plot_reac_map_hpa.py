# Standard library imports
from os import listdir
from os.path import isfile, join
import time
import pkg_resources as pkgr
#Third-party imports
import numpy as np
from itertools import product
# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.graph_module import *
from pomsimulator.modules.plotting_module import *
from pomsimulator.modules.helper_module import generate_CRN,generate_graphs
from configparser import ConfigParser


def main():
    Print_logo()
    ######################### User parameters ##############################################
    config_file = pkgr.resource_filename(__name__, "../inputs/config_PMo.pomsim")
    config = ConfigParser()
    config.read(config_file)
    # Input/output files
    system = config["Preparation"]["POM_system"]
    ADF_folder = pkgr.resource_filename(__name__, config["Preparation"]["adf_inputs_dir"])
    output_path = pkgr.resource_filename(__name__, config["Preparation"]["output_path"])

    isomorphism_matrix = output_path + "/np_IM_%s.csv" % system
    # Chemical parameters
    use_isomorphisms = config["Simulation"].getboolean("use_isomorphism")
    energy_threshold = float(config["Simulation"]["energy_threshold"])  # Maximum value for reaction energies
    proton_numb = int(config["Simulation"]["proton_numb"])  # Maximum difference in proton number of to species to react
    reference = config["Simulation"]["reference_types"].split(",")  # Reaction types of the network
    I, C0, temp = [float(config["Simulation"][prop]) for prop in ["I", "C0", "temp"]]
    min_pH, max_pH, step_pH = [float(config["Simulation"][prop]) for prop in ["min_pH", "max_pH", "step_pH"]]
    ref_compounds = config["Simulation"]["ref_compound"].split(",")
    # Internal parameters: These parameters are not meant to be routinely modified
    internal_conditions = config["InternalConditions"]
    internal_conditions["proton_numb"] = config["Simulation"]["proton_numb"]
    # CRN parameters
    full_CRN = config["CRN"].getboolean("Full_CRN")
    selected_model = int(config["CRN"]["Selected_model"])
    plot_3d = config["CRN"].getboolean("Plot_3D")

    # 1) Get ADF outputs

    print("1) Get ADF outputs")

    adf_files = sorted([ADF_folder + "/" + f for f in listdir(ADF_folder) if isfile(join(ADF_folder, f))])

    # 2) Graph creation: node=atom edge=bond
    print("2) Graph creation: node=atom edge=bond")

    G1_list, G1_labels, graphs_info = generate_graphs(adf_files, ref_compounds, system)
    stoich = [Lab_to_stoich(lab) for lab in G1_labels]

    # 3) Isomorphism matrix generation
    print("3) Isomorphism matrix generation")

    if use_isomorphisms:  # Calculate all isomorphisms
        diagonal = read_diagonal(isomorphism_matrix)
    else:  # Assume that all species are isomorphic
        diagonal = np.tri(graphs_info["num_molec"], graphs_info["num_molec"], 0)

    reac_idx, reac_energy, reac_type = Isomorphism_to_ChemicalReactions(G1_list, diagonal, graphs_info["water"], reference, system,
                                                                        energy_threshold, internal_conditions)

    print("3.1) NUMBER OF REACTIONS:", [len(list(map(len, reac_idx[i]))) for i in range(len(reac_idx))])
    R_idx, R_ene, R_type = sort_by_type(reac_idx, reac_energy, reac_type, graphs_info["compounds_set"], graphs_info["unique_labels"],
                                        G1_labels, exclude=ref_compounds)

    # 4) Speciation
    print("4) Generating CRN")


    plot_dict_details = dict(node_color='black', x_axis_lab='Metal/oxygen ratio', y_axis_lab='Number of protons', z_axis_lab='Number of metals',
                     plot_title='Reaction Map', colormap='RdYlGn_r', figsize=(8,6))

    number_models = np.prod([len(item) for item in R_type])
    print("4.1) Total number of models: ", number_models)

    plotting_dict = dict(full=full_CRN,dimension_3d=plot_3d,mod_idx=selected_model)
    crn_obj = generate_CRN(G1_list,stoich,reac_idx,reac_energy,reac_type,R_idx,R_ene,R_type,
                           plotting_dict,plot_dict_details)
    plt.savefig(output_path + "/Reaction_map_%s.png"%system,dpi=300)
    plt.show()
    print("Normal termination")
if __name__ == '__main__':
    main()
