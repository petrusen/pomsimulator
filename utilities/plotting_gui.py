# Standard library imports
import os
import re
from configparser import ConfigParser
import pkg_resources as pkgr
# Third-party imports
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for GUI threading
import matplotlib.pyplot as plt
# Local imports
from pomsimulator.modules.plotting_module import plot_speciation, get_color_phase_diagram
from pomsimulator.modules.helper_module import load_array, get_C0, phase_diagram_IPA, phase_diagram_HPA, get_config_to_dict
from pomsimulator.modules.text_module import Print_logo, Lab_to_Formula, Lab_to_stoich
from pomsimulator.modules.DataBase import Col_Dict_PMo


def plot_speciation_run(config_dict):
    """
    Generate speciation diagram plot for both IPA and HPA systems.
    
    Args:
        config_dict (dict or str): Configuration dictionary containing parameters for the function.
    Returns:
        str: Message indicating the completion of the function.
    """
    try:
        Print_logo()
        run_by_gui = True
        
        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False
        
        # Extract parameters
        raw_output_path = config_dict["Preparation"]["output_path"]
        
        # Handle path resolution more robustly
        if os.path.isabs(raw_output_path):
            # If it's already an absolute path, use it as-is
            output_path = raw_output_path
        else:
            # For relative paths, resolve relative to the utilities directory
            # Get the directory containing this module (utilities)
            current_dir = os.path.dirname(os.path.abspath(__file__))
            # Resolve the relative path from utilities directory
            output_path = os.path.join(current_dir, raw_output_path)
            output_path = os.path.normpath(output_path)
        
        system = config_dict["Preparation"]["POM_system"]
        
        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""
        
        # Get NPZ file path
        npz_file = config_dict["Speciation"].get("npz_file", "")
        print(npz_file)
        if not npz_file:
            npz_file = os.path.join(output_path, f"Array_{system}.npz")

        else:
            npz_file = output_path + "/" + config_dict["Speciation"]["npz_file"].split("/")[-1]
        # Normalize the path to resolve any .. components and duplications
        print("CCC", npz_file)
        # npz_file = os.path.normpath(npz_file)

        # Get m_idx (metal index for plotting)
        m_idx = int(config_dict["Speciation"].get("m_idx", "0"))
        
        # Get output image path
        output_img = os.path.join(output_path, f"Speciation_Diagram_{system}.png")
        
        # Check if we should stop
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"
        
        print(f"Loading speciation data from: {npz_file}")
        
        # Load the speciation array
        conc_arr, index_arr, C_ref, pH, speciation_labels = load_array(npz_file)
        C0 = get_C0(C_ref, m_idx=m_idx)
        conc_means = np.mean(conc_arr, axis=2)

        # Check if we should stop
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"
        
        # Get plot list (species to plot)
        plot_list = None
        if config_dict["Visualization"].get("plot_list", "all") != "all":
            plot_list = config_dict["Visualization"]["plot_list"].split(",")
            plot_list = [label.strip() for label in plot_list]

        # Get color dictionary
        col_dict = None
        col_dict_name = config_dict["Visualization"].get("col_dict", "")
        if col_dict_name and col_dict_name != "None (Default)":
            # Import color dictionaries from database
            try:
                from pomsimulator.modules.DataBase import color_dictionaries
                col_dict = color_dictionaries.get(col_dict_name)
            except (ImportError, AttributeError):
                print(f"Warning: Could not load color dictionary '{col_dict_name}'")

        print("Generating speciation diagram...")
        
        # Create the plot
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 4), constrained_layout=True)
        plot_speciation(conc_means, speciation_labels, pH, C0, plot_list, ax=ax, m_idx=m_idx, col_dict=col_dict)

        # Format legend with chemical formulas
        handle, labels = ax.get_legend_handles_labels()
        print(handle,labels)
        ax.legend(handle, [Lab_to_Formula(lab) for lab in labels])

        # Check if we should stop before saving
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            plt.close(fig)
            return "Stopped"

        # Save the plot
        plt.savefig(output_img, dpi=300)
        print(f"Speciation diagram saved to: {output_img}")
        # Only show plot if not running from GUI
        if not run_by_gui:
            plt.show()
        else:

            plt.close(fig)
        
        return "Normal termination"
        
    except Exception as e:
        print(f"An error occurred during speciation diagram generation: {str(e)}")
        return "Error"


def plot_phase_ipa_run(config_dict):
    """
    Generate IPA phase diagram plot.
    
    Args:
        config_dict (dict or str): Configuration dictionary containing parameters for the function.
    Returns:
        str: Message indicating the completion of the function.
    """
    try:
        Print_logo()
        run_by_gui = True
        
        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False
        
        # Extract parameters
        raw_output_path = config_dict["Preparation"]["output_path"]
        
        # Handle path resolution more robustly
        if os.path.isabs(raw_output_path):
            # If it's already an absolute path, use it as-is
            output_path = raw_output_path
        else:
            # For relative paths, resolve relative to the utilities directory
            # Get the directory containing this module (utilities)
            current_dir = os.path.dirname(os.path.abspath(__file__))
            # Resolve the relative path from utilities directory
            output_path = os.path.join(current_dir, raw_output_path)
            output_path = os.path.normpath(output_path)
        
        system = config_dict["Preparation"]["POM_system"]
        
        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""
        
        # Get phase directory or npz_info file
        phase_dir = config_dict["Speciation"].get("phase_dir", f"phase_diagram_{system}")
        npz_info_file = os.path.join(output_path, f"npz_info_{system}.dat")
        
        # Get species labels
        speciation_labels = config_dict["Speciation"].get("speciation_labels", "all").split(",")
        labels_file = os.path.join(output_path, f"labels_{system}.txt")
        if speciation_labels[0] == "all":
            try:
                with open(labels_file, "r") as flab:
                    speciation_labels = [item.strip() for item in flab.readlines()]
            except FileNotFoundError:
                print(f"Warning: Labels file not found: {labels_file}")
                return "Error"
        
        # Get plot list (species to plot)
        plot_list = None
        if config_dict["Visualization"].get("plot_list", "all") != "all":
            plot_list = config_dict["Visualization"]["plot_list"].split(",")
            plot_list = [label.strip() for label in plot_list]
        
        # Get color dictionary
        col_dict = None
        col_dict_name = config_dict["Visualization"].get("col_dict", "")
        if col_dict_name and col_dict_name != "None (Default)":
            try:
                from pomsimulator.modules.DataBase import color_dictionaries
                col_dict = color_dictionaries.get(col_dict_name)
            except (ImportError, AttributeError):
                print(f"Warning: Could not load color dictionary '{col_dict_name}'")
        
        # Check if we should stop
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"
        
        print(f"Loading phase diagram data from: {npz_info_file}")
        
        # Load the array with all speciations
        try:
            with open(npz_info_file, "r") as infile:
                npz_paths = [line.strip() for line in infile.readlines()]
        except FileNotFoundError:
            print(f"Error: NPZ info file not found: {npz_info_file}")
            return "Error"
        
        # Convert labels to stoichiometry for phase diagram calculation
        v_ctt = [Lab_to_stoich(lab) for lab in speciation_labels]
        
        # Generate phase diagram
        phase_diagram, C_list, pH = phase_diagram_IPA(npz_paths, v_ctt)
        
        # Check if we should stop
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"
        
        logC_list = np.log10(C_list)
        
        print("Generating IPA phase diagram...")
        
        # Get colors for phase diagram
        color_array, legend_elements = get_color_phase_diagram(phase_diagram, speciation_labels, col_dict)
        
        # Build the plot
        fig = plt.figure(constrained_layout=True, figsize=(6, 4))
        ax = fig.add_subplot()
        obj = ax.imshow(color_array,
                       extent=[max(pH), min(pH), min(logC_list), max(logC_list)],
                       origin='lower',
                       aspect='auto', alpha=1, interpolation=None)
        
        ax.set_xlabel('pH')
        ax.set_ylabel(r'$\log_{10}$' + ' Concentration (Molar)')
        ax.set_xlim(min(pH), max(pH))
        leg = ax.legend(handles=legend_elements, loc="center")
        leg.set_draggable(True)
        
        # Check if we should stop before saving
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            plt.close(fig)
            return "Stopped"
        
        # Path for saving the img
        output_img_path = re.sub("array.*npz", "", npz_paths[0]) + f"phase_diagram_{system}.png"
        plt.savefig(output_img_path, dpi=300)
        print(f"IPA phase diagram saved to: {output_img_path}")
        
        # Only show plot if not running from GUI
        if not run_by_gui:
            plt.show()
        else:
            plt.close(fig)
        
        return "Normal termination"
        
    except Exception as e:
        print(f"An error occurred during IPA phase diagram generation: {str(e)}")
        return "Error"


def plot_phase_hpa_run(config_dict):
    """
    Generate HPA phase diagram plot.
    
    Args:
        config_dict (dict or str): Configuration dictionary containing parameters for the function.
    Returns:
        str: Message indicating the completion of the function.
    """
    try:
        Print_logo()
        run_by_gui = True
        
        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False
        
        # Extract parameters
        raw_output_path = config_dict["Preparation"]["output_path"]
        
        # Handle path resolution more robustly
        if os.path.isabs(raw_output_path):
            # If it's already an absolute path, use it as-is
            output_path = raw_output_path
        else:
            # For relative paths, resolve relative to the utilities directory
            # Get the directory containing this module (utilities)
            current_dir = os.path.dirname(os.path.abspath(__file__))
            # Resolve the relative path from utilities directory
            output_path = os.path.join(current_dir, raw_output_path)
            output_path = os.path.normpath(output_path)
        
        system = config_dict["Preparation"]["POM_system"]
        
        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""
        
        # Get phase directory or npz_info file
        phase_dir = config_dict["Speciation"].get("phase_dir", f"phase_diagram_{system}")
        npz_info_file = os.path.join(output_path, f"npz_info_{system}.dat")
        
        # Get species labels
        speciation_labels = config_dict["Speciation"].get("speciation_labels", "all").split(",")
        labels_file = os.path.join(output_path, f"labels_{system}.txt")
        if speciation_labels[0] == "all":
            try:
                with open(labels_file, "r") as flab:
                    speciation_labels = [item.strip() for item in flab.readlines()]
            except FileNotFoundError:
                print(f"Warning: Labels file not found: {labels_file}")
                return "Error"
        
        # Get plot list (species to plot)
        plot_list = None
        if config_dict["Visualization"].get("plot_list", "all") != "all":
            plot_list = config_dict["Visualization"]["plot_list"].split(",")
            plot_list = [label.strip() for label in plot_list]
        
        # Get color dictionary - default to Col_Dict_PMo for HPA
        col_dict = Col_Dict_PMo  # Default for HPA
        col_dict_name = config_dict["Visualization"].get("col_dict", "")
        if col_dict_name and col_dict_name != "None (Default)":
            try:
                from pomsimulator.modules.DataBase import color_dictionaries
                user_col_dict = color_dictionaries.get(col_dict_name)
                if user_col_dict:
                    col_dict = user_col_dict
            except (ImportError, AttributeError):
                print(f"Warning: Could not load color dictionary '{col_dict_name}', using Col_Dict_PMo")
        
        # Check if we should stop
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"
        
        print(f"Loading phase diagram data from: {npz_info_file}")
        
        # Load the array with all speciations
        try:
            with open(npz_info_file, "r") as infile:
                npz_paths = [line.strip() for line in infile.readlines()]
        except FileNotFoundError:
            print(f"Error: NPZ info file not found: {npz_info_file}")
            return "Error"
        
        # Convert labels to stoichiometry for phase diagram calculation
        v_ctt = [Lab_to_stoich(lab) for lab in speciation_labels]
        
        # Generate phase diagram
        phase_diagram_X, phase_diagram_M, Ratio_list, pH = phase_diagram_HPA(npz_paths, v_ctt)
        
        # Check if we should stop
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"
        
        print("Generating HPA phase diagram...")
        
        # Build the plot with two subplots
        fig, ax = plt.subplots(ncols=2, nrows=1, constrained_layout=True, figsize=(8, 4))
        
        for pp, phase_diagram in enumerate([phase_diagram_X, phase_diagram_M]):
            color_array, legend_elements = get_color_phase_diagram(phase_diagram, speciation_labels, col_dict)
            obj = ax[pp].imshow(color_array,
                               extent=[max(pH), min(pH), min(Ratio_list), max(Ratio_list)],
                               origin='lower',
                               aspect='auto', alpha=1, interpolation=None)
            
            ax[pp].set_xlabel('pH')
            ax[pp].set_ylabel('M/X Ratio')
            ax[pp].set_xlim(min(pH), max(pH))
            leg = ax[pp].legend(handles=legend_elements, loc="center")
            leg.set_draggable(True)
        
        # Check if we should stop before saving
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            plt.close(fig)
            return "Stopped"
        
        # Path for saving the img
        output_img_path = re.sub("array.*npz", "", npz_paths[0]) + f"phase_diagram_{system}.png"
        plt.savefig(output_img_path, dpi=300)
        print(f"HPA phase diagram saved to: {output_img_path}")
        
        # Only show plot if not running from GUI
        if not run_by_gui:
            plt.show()
        else:
            plt.close(fig)
        
        return "Normal termination"
        
    except Exception as e:
        print(f"An error occurred during HPA phase diagram generation: {str(e)}")
        return "Error"