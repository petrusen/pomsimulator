# Standard library imports
from os import listdir
from os.path import isfile, join

# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.DataBase import *
from configparser import ConfigParser

def main():
    config_file = "../inputs/config_W.pomsim"
    config = ConfigParser()
    config.read(config_file)

    #### DEFINE VARIABLES ####

    ADF_folder = config["Preparation"]["adf_inputs_dir"]
    mol_folder = config["Preparation"]["mol_folder"]


    # 1) Get ADF outputs
    Print_logo()
    print("1) Get ADF outputs")

    adf_files = sorted([ADF_folder + f for f in listdir(ADF_folder) if isfile(join(ADF_folder, f))])

    # 2) Mol-File Generation:
    print("2) Mol-File Generation:")

    G1_list, water, G1_labels = list(), dict(), list()
    for idx, f in enumerate(adf_files):
        adf_dict = Bader_Parser(f)
        label = adf_dict['label']
        if label in ['H3O', 'H2O', 'H5O2', 'H4O2']:
            water[label] = adf_dict['Gibbs']
        else:
            write_molfile(mol_folder,limit_bonds=True, **adf_dict) # Creates .mol files from Bader connectivity

if __name__ == '__main__':
    main()
