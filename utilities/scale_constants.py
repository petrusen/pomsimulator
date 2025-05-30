# Standard library imports
from configparser import ConfigParser
import pkg_resources as pkgr

# Local imports
from pomsimulator.modules.text_module import Print_logo
from pomsimulator.modules.helper_module import LinearScaling
from pomsimulator.modules.stats_module import get_boxplot_data
from pomsimulator.modules.DataBase import *


def main():
    config_file = pkgr.resource_filename(__name__, "../inputs/config_W.pomsim")
    config = ConfigParser()
    config.read(config_file)

    system = config["Preparation"]["POM_system"]
    exp_set_name = config["Scaling"]["experimental_set"]
    ExpDict = experimental_constants[exp_set_name]
    output_path = pkgr.resource_filename(__name__, config["Preparation"]["output_path"])

    lgkf_file = output_path + "/logkf_%s.csv" % system
    scaling_params_file = output_path + "/regression_output.csv"
    scaling_mode = config["Scaling"]["scaling_mode"]
    Print_logo()

    scaling_params_dict = LinearScaling(lgkf_file, ExpDict, scaling_mode=scaling_mode,
                  output_scaling=scaling_params_file,Metal=system,output_path=output_path)
    print(scaling_params_dict)
    print("Normal termination")
if __name__ == '__main__':
    main()
