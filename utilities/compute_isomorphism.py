# Standard library imports
from os import listdir,makedirs
from os.path import isfile, join
from configparser import ConfigParser
import pkg_resources as pkgr
#Third-party imports
import numpy as np
# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.graph_module import *


def main():

    config_file = pkgr.resource_filename(__name__, "../inputs/config_W.pomsim")
    config = ConfigParser()
    config.read(config_file)
    # 0) Define input variables #### ###################################################################################

    ### Paths
    system = config["Preparation"]["POM_system"]
    mol_folder = pkgr.resource_filename(__name__, config["Preparation"]["mol_folder"])
    output_path = pkgr.resource_filename(__name__, config["Preparation"]["output_path"])
    output_file =  output_path + "/np_IM_%s.csv"% system
    Cores = int(config["Isomorphism"]["cores"]) #number of cores for the subgraph isomorphism calculation (caution!)
    os.makedirs(output_path, exist_ok=True)
    ### Variables for the reaction network


    # 1) Get ADF outputs ############################################################################################

    Print_logo()
    print("1) Get ADF outputs and generate parameters' output", "".join(["=" for _ in range(100)]))
    mol_files = sorted([mol_folder + "/" +  f for f in listdir(mol_folder) if isfile(join(mol_folder, f))])


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
            G1_list.append(Gi)

    num_molec = len(G1_list)

    # 3) Isomorphism matrix generation ##############################################################################

    print("3) Isomorphism matrix", "".join(["=" for _ in range(100)]))

    diagonal = Molecular_Graphs_to_Isomorphic_Matrix(G1_list, np.tri(num_molec, num_molec, 0), cores=Cores)

    np.savetxt(output_file,diagonal, fmt = '%d',delimiter = ',')

if __name__ == '__main__':
    main()
