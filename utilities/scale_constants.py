# Standard library imports
from configparser import ConfigParser
import pkg_resources as pkgr

# Local imports
from pomsimulator.modules.text_module import Print_logo
from pomsimulator.modules.helper_module import LinearScaling
from pomsimulator.modules.stats_module import get_boxplot_data
from pomsimulator.modules.DataBase import *


def scale_constants(config_dict):

    system = config_dict["Preparation"]["POM_system"]
    exp_set_name = config_dict["Scaling"]["experimental_set"]
    ExpDict = experimental_constants[exp_set_name]
    output_path = pkgr.resource_filename(__name__, config_dict["Preparation"]["output_path"])

    lgkf_file = output_path + "/logkf_%s.csv" % system
    scaling_params_file = output_path + "/regression_output.csv"
    scaling_mode = config_dict["Scaling"]["scaling_mode"]
    Print_logo()

    scaling_params_dict = LinearScaling(lgkf_file, ExpDict, scaling_mode=scaling_mode,
                  output_scaling=scaling_params_file,Metal=system,output_path=output_path)

    print(scaling_params_dict)
    print("Normal termination")

    return scaling_params_dict
