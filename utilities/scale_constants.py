# Local imports
from pomsimulator.modules.text_module import Print_logo
from pomsimulator.modules.helper_module import LinearScaling
from pomsimulator.modules.DataBase import *


def main():
    metal = "W"
    Labels = Labels_W
    ExpDict = None
    output_path = "../outputs/W_data"

    lgkf_file = output_path + "/logkf_W.txt"
    scaling_params_file = output_path + "/regression_output.csv"
    scaling_mode = "universal"
    Print_logo()
    scaling_params_dict = LinearScaling(lgkf_file, Labels, ExpDict, scaling_mode=scaling_mode,
                  output_scaling=scaling_params_file,Metal=metal,output_path=output_path)
    print(scaling_params_dict)

if __name__ == '__main__':
    main()
