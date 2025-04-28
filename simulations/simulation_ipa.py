# Standard library imports
from os import listdir
from os.path import isfile, join
import time
import random
import pkg_resources as pkgr
# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.helper_module import *
from configparser import ConfigParser

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

def main():
    Print_logo()
    ######################### User parameters ##############################################
    config_file = pkgr.resource_filename(__name__, "../inputs/config_W.pomsim")
    config = ConfigParser()
    config.read(config_file)
    # Input/output files

    system = config["Preparation"]["POM_system"]
    ADF_folder = pkgr.resource_filename(__name__, config["Preparation"]["adf_inputs_dir"])
    mol_folder = pkgr.resource_filename(__name__, config["Preparation"]["mol_folder"])
    output_path = pkgr.resource_filename(__name__, config["Preparation"]["output_path"])

    isomorphism_matrix = pkgr.resource_filename(__name__, output_path + "/np_IM.csv")
    formation_constants_file = output_path + "/logkf_%s.csv" % system
    CRN_file = output_path + "/CRN_%s.txt" % system
    labels_file = output_path + "/labels_%s.txt" % system
    simulation_file = output_path + "/simulation_parameters_%s.txt" % system

    # Operation parameters -> all from config file
    cores = int(config["Simulation"]["cores"])
    batch_size = int(config["Simulation"]["batch_size"])
    sample_perc = float(config["Simulation"]["sample_perc"])
    sample_type = config["Simulation"]["sample_type"] # random / all
    # Chemical parameters
    use_isomorphisms = config["Simulation"].getboolean("use_isomorphism")
    energy_threshold = float(config["Simulation"]["energy_threshold"])  # Maximum value for reaction energies
    proton_numb = int(config["Simulation"]["proton_numb"])     # Maximum difference in proton number of to species to react
    reference = config["Simulation"]["reference_types"].split(",") # Reaction types of the network
    I, C0, temp = [float(config["Simulation"][prop]) for prop in ["I","C0","temp"]]
    min_pH, max_pH, step_pH = [float(config["Simulation"][prop]) for prop in ["min_pH","max_pH","step_pH"]]
    ref_compound = config["Simulation"]["ref_compound"]

    # Internal parameters: These parameters are not meant to be routinely modified

    internal_conditions = config["InternalConditions"]
    internal_conditions["proton_numb"] = config["Simulation"]["proton_numb"]

    # 1) Get ADF outputs ############################################################################################

    print("1) Get ADF outputs")
    adf_files = sorted([ADF_folder + "/" +  f for f in listdir(ADF_folder) if isfile(join(ADF_folder, f))])

    # 2) Graph creation #############################################################################################

    print("2) Graph creation: node=atom edge=bond")

    G1_list, G1_labels, graphs_info = generate_graphs(adf_files, ref_compound, system)
    with open(labels_file,"w") as flab:
        flab.write("\n".join(G1_labels))
    # 3) Isomorphism matrix generation ##############################################################################

    print("3) Isomorphism matrix and reactions generation")

    if use_isomorphisms:  # Calculate all isomorphisms
        diagonal = read_diagonal(isomorphism_matrix)
    else:  # Assume that all species are isomorphic
        diagonal = np.tri(graphs_info["num_molec"], graphs_info["num_molec"], 0)

    reac_idx, reac_energy, reac_type = Isomorphism_to_ChemicalReactions(G1_list, diagonal, graphs_info["water"], reference, system,
                                                                        energy_threshold, internal_conditions)
    print("3.1) NUMBER OF REACTIONS:", [len(list(map(len, reac_idx[i]))) for i in range(len(reac_idx))])
    R_idx, R_ene, R_type = sort_by_type(reac_idx, reac_energy, reac_type, graphs_info["compounds_set"], graphs_info["unique_labels"],
                                        G1_labels, ref_compound)

    Write_Reactions(CRN_file, G1_labels, reac_idx, reac_type, reac_energy, stringreac_dict, molecularity_dict)

    # 4) Speciation #################################################################################################

    print("4) Solving MSCE Models and generate parameters' output")

    e_ctt = list(reac_energy[0])
    idx_ctt = list(reac_idx[0])
    type_ctt = list(reac_type[0])

    lgkf_params = dict(idx_ctt=idx_ctt, e_ctt=e_ctt, type_ctt=type_ctt, pH_grid=np.arange(min_pH, max_pH, step_pH),
                  init_guess=np.zeros(graphs_info["num_molec"]), I=I, C=C0, threshold=0.1, temp=temp)
    lgkf_params.update({k:graphs_info[k] for k in ["z_ctt","v_ctt","ref_idx"]})

    number_models = np.prod([len(item) for item in R_type])

    ### Printing output
    kwargs_input = dict()
    obj_list = [ADF_folder, mol_folder, formation_constants_file, CRN_file, simulation_file, cores, use_isomorphisms,
                energy_threshold, proton_numb, reference, I, C0, (min_pH, max_pH), step_pH, number_models, ref_compound,G1_labels]
    for s, o in zip(simulation_parameters_strings, obj_list):
        kwargs_input[s] = o
    write_simulationparameters(kwargs_input)

    print("4.1) Total number of models: ", number_models)

    start_time = time.time()

    mod_idx_vals = models_sampling(sample_type,number_models,sample_perc=sample_perc)
    data = compute_lgkf_loop(R_idx,R_ene,R_type,mod_idx_vals,number_models,lgkf_params,batch_size=batch_size,cores=cores)

    # 5) Writing Output #############################################################################################

    print("5) Create Output File")

    with open(formation_constants_file, "w") as out:
        header_str = "mod_idx," + ",".join(G1_labels)+"\n"
        out.write(header_str)
        for r,d in zip(mod_idx_vals,data):
            out.write(str(r)+","+",".join([str(di) for di in d])+'\n')
    print("Normal Termination. Execution time: " + str(round((time.time() - start_time), 4)) + " sec.")

if __name__ == '__main__':
    main()

