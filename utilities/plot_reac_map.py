# Standard library imports
from os import listdir
from os.path import isfile, join
import numpy as np
from itertools import product
import time

# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.graph_module import *
from pomsimulator.modules.plotting_module import *

def main():
    #### DEFINE VARIABLES ####

    ADF_Folder =  "../inputs/W_Set_PBE/"
    isomorphism_matrix =  "../utilities/np_IM.csv"

    target_model = 160
    ######################    CHEMICAL PARAMETERS    ############################
    use_isomorphisms = True
    energy_threshold = 30  # Maximum value for reaction energies
    proton_numb = 0  # Maximum difference in proton number of to species to react
    reference = ["P", "H2Ow1", 'H2Ow2', 'Cw1', 'Cw2', 'Cw3', 'Cw4', "A", 'HO', 'H3O']
    # isomorphism application
    ##############################################################################
    ######################    INTERNAL PARAMETERS   ##############################

    """ These parameters are not meant to be routinely modified """

    conditions_dict = {"proton_numb":proton_numb,"restrain_addition":2,"restrain_condensation":2,
                       "include_dimerization":True,"force_stoich":[11],"adjust_protons_hydration":True}

    # 1) Get ADF outputs
    Print_logo()
    print("1) Get ADF outputs")

    adf_files = sorted([ADF_Folder + f for f in listdir(ADF_Folder) if isfile(join(ADF_Folder, f))])

    # 2) Graph creation: node=atom edge=bond
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
            # write_molfile(MOL_Folder, Z_dict_inv, **adf_dict) # Creates .mol files from Bader connectivity

    print("Length:", len(G1_labels), G1_labels.index('W01O04-0H'), G1_labels)
    Z, valence, element = Determine_IPA(adf_dict['Z'])
    stoich = Molecule_Stoichiometry(G1_list, Z)
    compounds_set, unique_labels = Create_Stoich(G1_labels)

    # 3) Isomorphism matrix generation
    print("3) Isomorphism matrix generation")
    num_molec = len(G1_list)
    if use_isomorphisms:  # Calculate all isomorphisms
        diagonal = read_diagonal(isomorphism_matrix)
    else:  # Assume that all species are isomorphic
        diagonal = np.tri(num_molec, num_molec, 0)


    Reac_idx, Reac_energy, Reac_type = Isomorphism_to_ChemicalReactions(G1_list, diagonal, water, reference, "W",
                                                                        energy_threshold, conditions_dict)

    # 4) Speciation
    print("4) Solving MSCE Models")
    acc = 0
    e_ctt = list(Reac_energy[0])
    idx_ctt = list(Reac_idx[0])
    type_ctt = list(Reac_type[0])
    print(stoich)
    exclude = 'W01O04-0H'

    plot_dict = dict(node_color='black', x_axis_lab='Metal/oxygen ratio', y_axis_lab='Number of metals', z_axis_lab='Z Axis',
                     plot_title='Reaction Map', colormap='RdYlGn_r', figsize=(8,6))


    R_idx, R_ene, R_type = sort_by_type(Reac_idx, Reac_energy, Reac_type, compounds_set, unique_labels, G1_labels,
                                        exclude)
    print("NUMBER OF REACTIONS:", [len(list(map(len, Reac_idx[i]))) for i in range(len(Reac_idx))])
    begin = time.time()
    number_models = 0
    for _ in product(*R_type):
        number_models += 1

    print("Total number of models: ", number_models)
    end = time.time()
    print("Timing", round(end - begin, 2))


    comb_idx, comb_e, comb_type = product(*R_idx), product(*R_ene), product(*R_type)
    for idx_var, e_var, type_var in zip(comb_idx, comb_e, comb_type):
        if acc == target_model:
            '''Idx_new, e_new and type new consist of the reactions of a given speciation model, including acid-base reactions.
            Otherwise, the totality of the reactions is inside Reac_idx,Reac_energy and Reac_type. If the flag all models is True,
            you need to pass all the reactions. On the other hand, if the flag is false, only the new parameters have to be given.
            Plotting parameters as colormap, axis titles and plot title can be modified in the plot_dict above'''
            idx_new = idx_ctt + list(idx_var)
            e_new = e_ctt + list(e_var)
            type_new = type_ctt + list(type_var)
            fig,ax = Reaction_Map_2D_monometal(G1_list, idx_new, e_new, type_new, stoich, All_models=False,
                                       ploting_details_dict=plot_dict)
            #fig,axdict = Reaction_Map_3D_monometal(G1_list, Reac_idx,Reac_energy,Reac_type, stoich, All_models=True,
            #                                       ploting_details_dict=plot_dict)
            plt.savefig("../outputs/Reac_map_test.png",dpi=300)
            plt.show()
            break
        acc += 1

if __name__ == '__main__':
    main()
