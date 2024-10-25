# Local imports
from pomsimulator.modules.text_module import Print_logo
from pomsimulator.modules.helper_module import LinearScaling
from pomsimulator.modules.DataBase import *
import pandas as pd


def main():
    metal = "PMo"
    labels = Labels_PMo
    ExpDict = Pettersson_3I
    output_path = "../outputs/PMo_data"

    lgkf_file = output_path + "/logkf_PMo.txt"
    scaling_params_file = output_path + "/regression_output.csv"
    scaling_mode = "average"
    Print_logo()

    scaling_params_dict = LinearScaling(lgkf_file, labels, ExpDict, scaling_mode=scaling_mode,
                  output_scaling=scaling_params_file,Metal=metal,output_path=output_path)
    print(scaling_params_dict)

if __name__ == '__main__':
    main()
