
# Standard library imports
from os import listdir
from os.path import isfile, join
import time
import random
# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.helper_module import *

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

def main():
    Print_logo()
    ######################### User parameters ##############################################
    # Input/output files
    ADF_folder = "../inputs/PMo_testing/"
    mol_folder = "../inputs/PMo_testing_molfiles/"
    isomorphism_matrix = "../utilities/PMo_np_IM.csv"
    output_path = "../outputs/PMo_data"

    formation_constants_file = output_path + "/logkf_PMo_test.txt"
    CRN_file = output_path + "/PMo_CRN.txt"
    simulation_file = output_path + "/simulation_parameters_PMo.txt"

    # Operation parameters
    cores = 20
    batch_size = 100
    sample_perc = 10
    sampl_type = "random" # random / all

    # Chemical parameters
    use_isomorphisms = False
    energy_threshold = 15  # Maximum value for reaction energies
    proton_numb = 1     # Maximum difference in proton number of to species to react
    reference = ["P", "H2Ow1", 'H2Ow2', 'Cw1', 'Cw2', 'Cw3', 'Cw4', "A", 'HO', 'H3O']  # Reaction types of the network
    I, C0 = 0.25, 0.005
    temp = 298.15
    min_pH, max_pH, grid = 0, 35, 70
    ref_compounds = ['P01Mo00O04-0H','P00Mo01O04-0H']
    POM_system = "P_Mo" # User should input the system as X_M (X for the heteroatom and M for the metal)
    # Internal parameters: These parameters are not meant to be routinely modified

    conditions_dict = {"proton_numb":proton_numb,"restrain_addition":12,"restrain_condensation":12,
                       "include_dimerization":True,"force_stoich":[],"adjust_protons_hydration":True}

    # 1) Get ADF outputs ############################################################################################

    print("1) Get ADF outputs")
    adf_files = sorted([ADF_folder + f for f in listdir(ADF_folder) if isfile(join(ADF_folder, f))])

    # 2) Graph creation #############################################################################################

    print("2) Graph creation: node=atom edge=bond")

    G1_list, G1_labels, graphs_info = generate_graphs(adf_files, ref_compounds,POM_system)

    # 3) Isomorphism matrix generation ##############################################################################

    print("3) Isomorphism matrix and reactions generation")

    if use_isomorphisms:  # Calculate all isomorphisms
        diagonal = read_diagonal(isomorphism_matrix)
    else:  # Assume that all species are isomorphic
        diagonal = np.tri(graphs_info["num_molec"], graphs_info["num_molec"], 0)

    reac_idx, reac_energy, reac_type = Isomorphism_to_ChemicalReactions(G1_list, diagonal, graphs_info["water"], reference, POM_system,
                                                                        energy_threshold, conditions_dict)
    print("3.1) NUMBER OF REACTIONS:", [len(list(map(len, reac_idx[i]))) for i in range(len(reac_idx))])
    R_idx, R_ene, R_type = sort_by_type(reac_idx, reac_energy, reac_type, graphs_info["compounds_set"], graphs_info["unique_labels"],
                                        G1_labels, exclude=ref_compounds)

    Write_Reactions(CRN_file, G1_labels, reac_idx, reac_type, reac_energy, stringreac_dict, molecularity_dict)

    # 4) Speciation #################################################################################################

    print("4) Solving MSCE Models and generate parameters' output")

    e_ctt = list(reac_energy[0])
    idx_ctt = list(reac_idx[0])
    type_ctt = list(reac_type[0])

    lgkf_params = dict(idx_ctt=idx_ctt, e_ctt=e_ctt, type_ctt=type_ctt,
                       pH_grid=np.linspace(min_pH, max_pH, grid),
                  init_guess=np.zeros(graphs_info["num_molec"]), I=I, C_X=C0,C_M=C0, threshold=0.1,temp=temp)
    lgkf_params.update({k:graphs_info[k] for k in ["z_ctt","v_ctt","ref_idx"]})

    begin = time.time()
    number_models = 0
    for _ in product(*R_type):
        number_models += 1
    end = time.time()

    ### Printing output
    kwargs_input = dict()
    obj_list = [ADF_folder, mol_folder, formation_constants_file, CRN_file, simulation_file, cores, use_isomorphisms,
                energy_threshold, proton_numb, reference, I, C0, (min_pH, max_pH), grid, number_models, ref_compounds,G1_labels]
    for s, o in zip(simulation_parameters_strings, obj_list):
        kwargs_input[s] = o
    write_simulationparameters(kwargs_input)

    print("4.1) Total number of models: ", number_models)
    print("Timing", round(end - begin, 2))

    start_time = time.time()

    mod_idx_vals = models_sampling(sampl_type,number_models,sample_perc=sample_perc)
    data = compute_lgkf_loop(R_idx,R_ene,R_type,mod_idx_vals,number_models,lgkf_params,batch_size=batch_size,cores=cores)


    # 5) Writing Output #############################################################################################

    print("5) Create Output File")

    with open(formation_constants_file, "w") as out:
        for r,d in zip(mod_idx_vals,data):
            out.write(str(r)+","+",".join([str(di) for di in d])+'\n')
    print("Normal Termination. Execution time: " + str(round((time.time() - start_time), 4)) + " sec.")

if __name__ == '__main__':
    main()

