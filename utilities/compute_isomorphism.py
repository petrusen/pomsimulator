# Standard library imports
from os import listdir
from os.path import isfile, join
import numpy as np

# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.graph_module import *


def main():
    # 0) Define input variables #### ###################################################################################

    ### Paths

    MOL_Folder =  "../inputs/W_Set_PBE_molfiles/"
    output_path = "../outputs/W_data"
    output_file =  output_path + "/np_IM_W.csv"
    Cores = 5  # cpu_count()  #number of cores for the simulation (caution!)

    ### Variables for the reaction network


    # 1) Get ADF outputs ############################################################################################

    Print_logo()
    print("1) Get ADF outputs and generate parameters' output", "".join(["=" for _ in range(100)]))
    mol_files = sorted([MOL_Folder + f for f in listdir(MOL_Folder) if isfile(join(MOL_Folder, f))])


    # 2) Graph creation #############################################################################################

    print("2) Graph creation", "".join(["=" for _ in range(100)]))
    G1_list = list()
    for idx, f in enumerate(mol_files):
        mol_dict = Mol_Parser_2(f)
        label = mol_dict['label']
        if label in ['H3O', 'H2O', 'H5O2', 'H4O2']:
            continue
        else:
            Gi = Molecule_to_Graph_from_molfile(idx, mol_dict["Z"], mol_dict["bonds"], mol_dict["label"])
            #print(Gi,Gi.edges)
            G1_list.append(Gi)

    num_molec = len(G1_list)

    # 3) Isomorphism matrix generation ##############################################################################

    print("3) Isomorphism matrix", "".join(["=" for _ in range(100)]))

    diagonal = Molecular_Graphs_to_Isomorphic_Matrix(G1_list, np.tri(num_molec, num_molec, 0), cores=Cores)

    np.savetxt(output_file,diagonal, fmt = '%d',delimiter = ',')

if __name__ == '__main__':
    main()
