# Standard library imports
from os import listdir
from os.path import isfile, join
import numpy as np
from multiprocessing import Pool, cpu_count
from itertools import repeat, product, islice, compress
import time

# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.graph_module import *
from pomsimulator.modules.msce_module import *
from pomsimulator.modules.DataBase import *


def main():

    #### DEFINE VARIABLES ####

    ADF_folder = "../inputs/PMo_testing/"
    mol_folder = "../inputs/W_Set_PBE_molfiles/"
    isomorphism_matrix = "../utilities/np_IM.csv"
    formation_constants_file = "../outputs/logkf_PMo.txt"
    CRN_file = "../outputs/PMo_CRN.txt"
    simulation_file = "../outputs/simulation_parameters_PMo.txt"
    cores = 20
    batch_size = 200



    ######################    CHEMICAL PARAMETERS    ############################
    use_isomorphisms = False
    energy_threshold = 15  # Maximum value for reaction energies
    proton_numb = 1     # Maximum difference in proton number of to species to react
    reference = ["P", "H2Ow1", 'H2Ow2', 'Cw1', 'Cw2', 'Cw3', 'Cw4', "A", 'HO', 'H3O']  # Reaction types of the network
    POM = "P_Mo"
    ##############################################################################
    ######################    INTERNAL PARAMETERS   ##############################

    """ These parameters are not meant to be routinely modified """

    conditions_dict = {"proton_numb":proton_numb,"restrain_addition":12,"restrain_condensation":12,
                       "include_dimerization":True,"force_stoich":[],"adjust_protons_hydration":True}

    ### Variables for the speciation
    I, C0 = 0.25, 0.005
    min_pH, max_pH, grid = 0, 35, 70
    temp = 298.15
    ref_compounds = ['P01Mo00O04-0H','P00Mo01O04-0H']

    # 1) Get ADF outputs ############################################################################################

    Print_logo()
    print("1) Get ADF outputs")
    
    adf_files = sorted([ADF_folder + f for f in listdir(ADF_folder) if isfile(join(ADF_folder, f))])

    # 2) Graph creation #############################################################################################

    print("2) Graph creation: node=atom edge=bond")
    
    G1_list, water, G1_labels = list(), dict(), list() 
    for idx, f in enumerate(adf_files):
        adf_dict = Bader_Parser(f)
        label = adf_dict['label']
        if label in ['H3O', 'H2O', 'H5O2', 'H4O2']:
            water[label] = adf_dict['Gibbs']
        else:
            Gi = Molecule_to_Graph(idx, **adf_dict)
            G1_list.append(Gi)
            G1_labels.append(label)

    # print("Length:", len(G1_labels), G1_labels.index('W01O04-0H'), G1_labels)
    elements = POM.split("_")
    Z = [Z_dict[elem] for elem in elements]
    valence = [valence_dict[elem] for elem in elements]

    charges = Molecule_Charge(G1_list, Z, valence)
    stoich = Molecule_Stoichiometry(G1_list, Z)
    compounds_set, unique_labels = Create_Stoich(G1_labels)
    num_molec = len(G1_list)

    # 3) Isomorphism matrix generation ##############################################################################

    print("3) Isomorphism matrix generation")

    if use_isomorphisms:  # Calculate all isomorphisms
        diagonal = read_diagonal(isomorphism_matrix)
    else:  # Assume that all species are isomorphic
        diagonal = np.tri(num_molec, num_molec, 0)


    reac_idx, reac_energy, reac_type = Isomorphism_to_ChemicalReactions(G1_list, diagonal, water, reference, POM,
                                                                        energy_threshold, conditions_dict)

    Write_Reactions(CRN_file, G1_labels, reac_idx, reac_type, reac_energy, stringreac_dict, molecularity_dict)

    # 4) Speciation #################################################################################################

    print("4) Solving MSCE Models and generate parameters' output")

    e_ctt = list(reac_energy[0])
    idx_ctt = list(reac_idx[0])
    type_ctt = list(reac_type[0])

    ref_idxs = [G1_labels.index(ref) for ref in ref_compounds]
    kwargs = dict(idx_ctt=idx_ctt, e_ctt=e_ctt, type_ctt=type_ctt, z_ctt=charges, v_ctt=stoich,
                  ref_idxs=ref_idxs, pH_grid=np.linspace(min_pH, max_pH, grid),
                  init_guess=np.zeros(num_molec), I=I, C_X=C0, C_M=C0, threshold=0.1, temp=temp)

    R_idx, R_ene, R_type = sort_by_type(reac_idx, reac_energy, reac_type, compounds_set, unique_labels, G1_labels,
                                        exclude=ref_compounds)

    print("4.1) NUMBER OF REACTIONS:", [len(list(map(len, reac_idx[i]))) for i in range(len(reac_idx))])

    begin = time.time()
    number_models = 0
    for _ in product(*R_type):
        number_models += 1

    print("4.2) Total number of models: ", number_models)
    end = time.time()
    print("Timing", round(end - begin, 2))

    data = list()
    start_time = time.time()

    low = 0
    up = number_models


    ### Printing output
    kwargs_input = dict()
    obj_list = [ADF_folder, mol_folder, formation_constants_file, CRN_file, simulation_file, cores, use_isomorphisms,
                energy_threshold, proton_numb, reference, I, C0, (min_pH, max_pH), grid, (low, up), ref_compounds, G1_labels]
    for s, o in zip(simulation_parameters_strings, obj_list):
        kwargs_input[s] = o
    write_simulationparameters(kwargs_input)

    mod_idx_vals = list(range(low, up))
    models_to_explore = set(mod_idx_vals)

    n_batches = int(len(mod_idx_vals) / batch_size)
    print("Number of batches = %d" % n_batches)

    _idx_var, _e_var, _type_var = product(*R_idx), product(*R_ene), product(*R_type)
    bool_sample = (idx in models_to_explore for idx in range(number_models))
    models_to_calculate = compress(zip(_idx_var, _e_var, _type_var),bool_sample)
    var = list()
    acc = 0
    for obj in models_to_calculate:
        var_args = (obj[0],obj[1],obj[2])
        var.append(var_args)
        acc += 1

    kwargs_iter = repeat(kwargs)

    print("4.3) Enter formation constant calculation")

    for idx in range(n_batches + 1):
        t0 = time.time()

        low_lim = batch_size * idx
        up_lim = batch_size * (idx + 1)
        islice_end = batch_size
        if idx == n_batches:
            current_var = var[low_lim:]
        else:
            current_var = var[low_lim:up_lim]
        args_iter = current_var

        with Pool(cores) as ThreadPool:  # HPC
            data = data + starmap_with_kwargs(ThreadPool, Speciation_from_Equilibrium_bimetal,
                                              args_iter, kwargs_iter)
        t1 = time.time()
        progress = idx * 100 / n_batches
        named_tuple = time.localtime()  # get struct_time
        time_string = time.strftime("%m/%d/%Y, %H:%M:%S", named_tuple)
        print(time_string,
              "[" + "".join(['#' if i < progress else " " for i in range(0, 100, 2)]) + "]" + " progress=%.3f" % progress,
              "time of batch = %.2f s" % (t1 - t0))

    # 5) Writing Output #############################################################################################

    print("5) Create Output File")

    with open(formation_constants_file, "w") as out:
        for r,d in enumerate(data):
            out.write(str(r)+","+",".join([str(di) for di in d])+'\n')
    print("Normal Termination. Execution time: " + str(round((time.time() - start_time), 4)) + " sec.")

if __name__ == '__main__':
    main()

