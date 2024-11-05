# Standard library imports
from os import listdir
from os.path import isfile, join

# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.DataBase import *


def main():
    #### DEFINE VARIABLES ####

    ADF_Folder = "../inputs/PMo_Set/"
    MOL_Folder = "../inputs/PMo_Set_molfiles/"


    # 1) Get ADF outputs
    Print_logo()
    print("1) Get ADF outputs")

    adf_files = sorted([ADF_Folder + f for f in listdir(ADF_Folder) if isfile(join(ADF_Folder, f))])

    # 2) Mol-File Generation:
    print("2) Mol-File Generation:")

    G1_list, water, G1_labels = list(), dict(), list()
    for idx, f in enumerate(adf_files):
        adf_dict = Bader_Parser(f)
        label = adf_dict['label']
        if label in ['H3O', 'H2O', 'H5O2', 'H4O2']:
            water[label] = adf_dict['Gibbs']
        else:
            write_molfile(MOL_Folder,limit_bonds=True, **adf_dict) # Creates .mol files from Bader connectivity

if __name__ == '__main__':
    main()
