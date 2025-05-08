# Standard library imports
import re
import pandas as pd
import datetime
import numpy as np

# Local imports
from pomsimulator.modules.DataBase import Z_dict,Z_dict_inv


def Lab_to_stoich(Label):
    '''From a molecular label of the form MxxOyy-zH, gathers the corresponding stoichiometry
     as a tuple of integers (xx,yy,z)
    Args:
        Label: string, name of the form MxxOyy-zH
    Returns:
        stoich: tuple of integers, (xx,yy,z)
    '''
    param = "([A-z]{1,2})([0-9]{1,2})"
    paramH = "([-])([0-9]{1,2})([H])"
    matches = re.findall(param,Label) + re.findall(paramH,Label)
    stoich = [int(item[1]) for item in matches]
    return stoich


def Lab_to_Formula(Label):
    '''From a molecular label of the form MxxOyy-zH, formats it to transform numbers to
    subindices and include the corresponding charges, assuming fully-oxidized central metals
    of groups 5 and 6. 
    Args:
        Label: string, name of the form MxxOyy-zH
    Returns:
        out_str: string, LaTeX-formatted formula for plotting
    '''
    charges = {"O": -2, "H": 1, "V": 5, "Nb": 5, "Ta": 5, "Mo": 6, "W": 6, "P": 5}
    param = "([A-z]{1,2})([0-9]{1,2})"
    paramH = "([-])([0-9]{1,2})([H])"
    # ([A-z]{1,2})([0-9]{1,2})([A-Z])([0-9]{1,2})([-][H])([0-9]{1,2})"
    metal_info = re.findall(param, Label)

    proton_info = re.findall(paramH, Label)
    elements = [item[2] for item in proton_info] + [item[0] for item in metal_info]
    coefs = [int(item[1]) for item in (proton_info + metal_info)]

    total_charge = sum([charges[elem] * cx for cx, elem in zip(coefs, elements)])

    out_str = r"$\mathrm{["
    for cx, elem in zip(coefs, elements):
        if cx == 0:
            continue
        elif cx == 1:
            out_str += elem
        else:
            out_str += "%s_{%d}" % (elem, cx)
    if total_charge == 1:
        out_str += "]^+}$"
    elif total_charge == -1:
        out_str += "]^-}$"
    elif total_charge < 0:
        out_str += "]^{%d\!-}}$" % abs(total_charge)
    elif total_charge > 0:
        out_str += "]^{%d\!+}}$" % abs(total_charge)
    else:
        out_str += "]}$"

    return out_str

def write_xyzfile(path, **kwargs):
    """
    Creates a cartesian file from the atomic numbers and the atom coordinates. It reduces the dependence
    of POMSimulator to external libraries.
    
    Args: 
        path: string, path to where the cartesian file will be written
        kwargs: dictionary, list of atomic numbers, cartesian coordinates, and file label

    Returns:
        bolean: True, if successful

    """

    Z, xyz, Label = kwargs['Z'], kwargs['xyz'], kwargs['label']
    with open(path + Label + ".xyz", "W_P") as outfile:
        outfile.write(str(len(Z)))
        outfile.write("\n")
        outfile.write("\n")
        for zi, xyz_i in zip(Z, xyz):
            outfile.write(Z_dict_inv[zi])
            outfile.write("  ")
            for xyz_i_j in xyz_i:
                outfile.write(str(xyz_i_j))
                outfile.write("   ")
            outfile.write("\n")
    
    return True

def write_molfile(path, verbose=True, limit_bonds=False, **kwargs):
    """
    Creates a mol file from the atomic numbers and the atom coordinates. It reduces the dependence
    of POMSimulator to external libraries.

    Args:
        path: string, path to where the cartesian file will be written
        verbose: bolean, prints out the progress
        limit_bonds: boolean, only allow atom - oxygen bonds (excluding O - O).
        kwargs: dictionary, list of atomic numbers, cartesian coordinates, chemical connectivity, and file label
    
    Returns
        bolean: True for succesfully created .mol files

    """

    Z, Bonds, xyz, Label = kwargs['Z'], kwargs['bonds'], kwargs['xyz'], kwargs['label']
    countline = " {a} {b}  0  0  0  0  0  0  0  0999 V2000\n"
    atomblock = "   {x}   {y}   {z} {Z}    0  0  0  0  0  0  0  0  0  0  0  0\n"
    bondblock = " {zi} {zj}  1  0  0  0  0\n"
   
    if len(Bonds) == 0: # if no chemical bonds, exit the function
        if verbose: print("WARNING: missing bond connectivity for {a}".format(a=Label));
        return False

    if limit_bonds:
        valid_bonds = []
        for zz in Bonds:
            zi, zj = [int(o) for o in zz]
            test = tuple(sorted([Z[zi-1],Z[zj-1]]))
            # only allow bonds involving one oxygen and one non-oxygen atom
            if test.count(8) == 1:
                valid_bonds.append(zz)
        Bonds = valid_bonds

    with open(path + '/' +Label + ".mol", 'w') as outfile:
        outfile.write(" \n") # Header - Title
        outfile.write(" Created with POMSimulator \n") # Header - Original Program
        outfile.write("\n") # Header - Black line
        strZ, strB = '{:>2.0f}'.format(len(Z)), '{:>2.0f}'.format(len(Bonds))
        outfile.write(countline.format(a=strZ, b=strB))
        for Z_i, xyz_i in zip(Z, xyz):
            xi, yi, zi = [float(o) for o in xyz_i]
            strx, stry, strz = '{:>7.4f}'.format(xi), '{:>7.4f}'.format(yi), '{:>7.4f}'.format(zi)
            outfile.write(atomblock.format(x=strx, y=stry, z=strz, Z=Z_dict_inv[float(Z_i)]))
        for zz in Bonds:
            zi, zj = [int(o) for o in zz]
            strzi, strzj =  '{:>2.0f}'.format(zi),  '{:>2.0f}'.format(zj)
            outfile.write(bondblock.format(zi=strzi, zj=strzj))
        outfile.write("M  END") # End   
        if verbose: print("Succesful generation of the mol file for {a}".format(a=Label))
        return True
        
def Mol_Parser(output_file):
    """
    Parser of .mol files which returns the cartesian coordinates, and the bond connectivity.
    
    Args:
        output_file: string, path to where the mol file is 
    """

    mol_dict = {}
    with open(output_file,"r") as out:
        acc = 0
        natoms, nbonds = 0, 0
        xyz, bonds = list(), list()
        for line in out:
            list_line = line.split()
            if acc == 3:
                natoms, nbonds = int(list_line[0]), int(list_line[1])
            elif acc > 3 and acc <= 3 + natoms:
                atom_i = list_line[3]
                xyz.append(Z_dict[atom_i])
            elif acc > 3 + natoms and acc <= 3 + natoms + nbonds:
                bond = (int(list_line[0]), int(list_line[1]))
                bonds.append(bond)
            acc += 1
    mol_dict['xyz'] = xyz
    mol_dict['bonds'] = bonds
    
    return mol_dict


def Mol_Parser_2(mol_file):
    """
    Parser of .mol files which returns the cartesian coordinates, and the bond connectivity.

    Args:
        mol_file: string, path to where the mol file is
    """

    mol_dict = {}
    with open(mol_file, "r") as out:
        label = mol_file.split("/")[-1].split(".")[0]
        acc = 0
        natoms, nbonds = 0, 0
        Z_list, bonds = list(), list()
        for line in out:
            list_line = line.split()
            if acc == 3:
                natoms, nbonds = int(list_line[0]), int(list_line[1])
            elif acc > 3 and acc <= 3 + natoms:
                atom_i = list_line[3]
                Z_list.append(Z_dict[atom_i])
            elif acc > 3 + natoms and acc <= 3 + natoms + nbonds:
                bond = (int(list_line[0]), int(list_line[1]))
                bonds.append(bond)
            acc += 1
    mol_dict['Z'] = Z_list
    mol_dict['bonds'] = bonds
    mol_dict['label'] = label



    return mol_dict

def Bader_Parser(output_file):
    """
    Parser of output files from Amsterdam Density Functional (ADF) software package [1].
    The calculation of corresponds to a geometry optimization, with frequencies and the
    added topological analysis of Bader [2].

    [1] te Velde, G.; Bickelhaupt, F. M.; Baerends, E. J.; Fonseca Guerra, C.;
        van Gisbergen, S. J. A.; Snijders, J. G.; Ziegler, T.
        Chemistry with ADF. J. Comput. Chem. 2001, 22 (9), 931–967.
    [2] Bader, R. F. W. A Quantum Theory of Molecular Structure and Its Applications.
        Chem. Rev. 1991, 91 (5), 893–928.

    Args:
        output_file: string, path to the output file of ADF2019.

    Returns:
        adf_dict: dictionary, cartesian coordinates, atomic numbers, chemical bonds, label, Gibbs and
        enthalpy energies, molecular charge, and bonding energy
    
    """

    adf_dict = {}
    with open(output_file,"r") as out:
        critical_points = dict([])
        connectivity = list()
        parsed_data = list()
        cartesian_coords, atomic_value = list(), list()
        pointer_conn = False
        pointer_xyz, pointer_xyz2 = False, 0
        pointer_BE, Bonding_Energy, Gibbs = True, 0, 0
        pointer_Z, second_pointer_Z, list_Z = False, False, list()
        pointer2_Z, second_pointer2_Z = False, False
        Enthalpy = 0
        for idx, line in enumerate(out):
            list_line = line.split()
            list_line2 = line.split(":")
            if "ATOMS" in list_line:
                pointer_Z = True
            elif "FRAGMENTS" in list_line:
                pointer_Z = False
            if pointer_Z == True:
                z_labels = Z_dict.keys()

                for zlab in z_labels:
                    if zlab in list_line:
                        if 'CHARGE' not in list_line:
                            list_Z.append(Z_dict[zlab])
                            x,y,z = list_line[2], list_line[3], list_line[4]
                            cartesian_coords.append([x,y,z])
            if "  THIS IS THE GEOMETRY IN THE CP SEARCH (ANGSTROM)" \
                    in list_line2 or pointer_xyz == True:  # Obtain xyz matrix and atom types
                pointer_xyz = True
                parsed_data.append(list_line)
            if "CRITICAL" in list_line:  # Obtain critical points
                pointer_xyz = False
                key = str()
                for i in range(len(list_line[0:5])):
                    key = key + list_line[i] + ' '
                critical_points[key] = int(list_line[5])
            if "#BP" in list_line or pointer_conn == True:  # Obtain connectivity
                pointer_conn = True
                if len(list_line) >= 2 and list_line[0] != '#BP':
                    connectivity.append((int(list_line[2]), int(list_line[3])))
                if len(list_line) == 0:
                    pointer_conn = False
            if 'Bonding' in list_line and pointer_BE == True:
                Bonding_Energy = float(list_line[5])
                pointer_BE = False
            if 'Enthalpy' in list_line and 'H:' in list_line:
                Enthalpy = float(list_line[4])
            if 'Gibbs' in list_line:
                Gibbs = float(list_line[5])
            if len(list_line) > 0 and 'charge' in list_line[0]:
                charge = int(float(list_line[1]))
            else:
                charge = 0

    labels = output_file.split("/")[-1].split(".")[0]
    adf_dict['xyz'] = cartesian_coords
    adf_dict['Z'] = list_Z
    adf_dict['bonds'] = connectivity
    adf_dict['label'] = labels
    adf_dict['Gibbs'] = Gibbs
    adf_dict['Enthalpy'] = Enthalpy
    adf_dict['charge'] = charge
    adf_dict['Electronic'] = Bonding_Energy
    
    return adf_dict


def Write_Reactions(path, G1_labels, all_reac_idx, all_reac_type, all_reac_e_eq, stringreac_dict, molecularity_dict):
    """
    Creates a text file with all the chemical reactions, and their respective reaction free energies.
    """

    with open(path, "w") as outfile:
        for reac_idx,reac_type,reac_e in zip(all_reac_idx, all_reac_type, all_reac_e_eq):
            for idx,rtype,ei in zip(reac_idx,reac_type,reac_e):
                e = str(round(ei,2))
                if molecularity_dict[rtype] == 1:
                    p, r1 = idx
                    R1, P = G1_labels[r1], G1_labels[p]
                    outfile.write(stringreac_dict[rtype].format(R1=G1_labels[r1], P=G1_labels[p], G=e))
                elif molecularity_dict[rtype] == 2:
                    p, r1, r2 = idx
                    R1, R2, P = G1_labels[r1], G1_labels[r2], G1_labels[p]
                    outfile.write(stringreac_dict[rtype].format(R1=R1, R2=R2, P=P, G=e))
    return True


def Read_csv(path):
    '''Helper function to read formation constant files, in CSV format, as a Pandas DataFrame
    already removing lines without valid constants and handling large files.
    Args:
        path: string, path to the file with formation constants
    Returns:
        my_data. Pandas.DataFrame with formation constants.

    '''
    my_data = pd.read_csv(path, sep=',', index_col=0, skip_blank_lines=True,
                          on_bad_lines='skip',na_values="None", low_memory=False)

    my_data.dropna(axis=0, how='all', inplace=True)

    #my_data.index -= 1
    return my_data

def Print_logo():
    banner = ''' 
    ---------------------------------------------------------------------------------------------------------------------------------
    |     *************************************************                                                                         |
    |     * ICIQ (Institut Català d'Investigació Química) *                                                                         |
    |     * Tarragona                                     *                                                                         |
    |     *************************************************                                                                         |
    |                                                                                                                               |
    |       8888888b.    .d88888b.   888b     d888  .d8888b.  d8b                        888          888                           |
    |       888   Y88b  d88P" "Y88b  8888b   d8888 d88P  Y88b Y8P                        888          888                           |
    |       888    888  888     888  88888b.d88888 Y88b.                                 888          888                           |
    |       888   d88P  888     888  888Y88888P888  "Y888b.   888 88888b.d88b.  888  888 888  8888b.  888888 .d88b.  888d888        |
    |       8888888P"   888     888  888 Y888P 888     "Y88b. 888 888 "888 "88b 888  888 888     "88b 888   d88""88b 888P"          |
    |       888         888     888  888  Y8P  888       "888 888 888  888  888 888  888 888 .d888888 888   838  888 888            |
    |       888         Y88b. .d88P  888   "   888 Y88b  d88P 888 888  888  888 Y88b 888 888 888  888 Y88b. Y88..88P 888            |
    |       888          "Y88888P"   888       888  "Y8888P"  888 888  888  888  "Y88888 888 "Y888888  "Y888 "Y88P"  888            |
    |                                                                                                                               |
    |                                                                 Created by  ENRIC PETRUS     2020                             |
    |                                                                 Authors: ENRIC PETRUS, MIREIA SEGADO-CENTELLAS, CARLES BO     |
    |                                                                 Developer team: ENRIC PETRUS, JORDI BUILS, DIEGO GARAY-RUIZ   |   
    ---------------------------------------------------------------------------------------------------------------------------------                                                                                                   
     '''
    print(banner)
    return None
                                                                                                                                                                                                                                                                                                                                                                  
def write_simulationparameters(kwargs):
    """
    Generates an output file with the variables and parameters chosen for a given simulation
    """

    with open(kwargs["Simulation Parameters File"], 'w') as outfile:


        outfile.write("=============================================================================================\n")
        outfile.write("·······································POMSimulator··········································\n")
        outfile.write("=============================================================================================\n")
        outfile.write("Information about the software's license, documentation and citation can be found at:   \n")
        outfile.write("https://gitlab.com/enricpp/pomsimulator_release                                         \n\n")
        outfile.write("Authors: Enric Petrus, Mireia Segado-Centellas, Carles Bo                                \n")
        outfile.write("Institute of Chemical Research of Catalonia (ICIQ) & Universitat Rovira i Virgili (URV) \n")
        outfile.write("Developer team: Enric Petrus, Jordi Buils, Diego Garay-Ruiz                                  \n")
        outfile.write("Contact Mail: enricpz@icloud.com                                                             \n")
        outfile.write("=============================================================================================\n\n")

        outfile.write("================\n")
        outfile.write("System Settings:\n")
        outfile.write("================\n")
        l_paths = ['ADF Folder', 'MOL Folder', "Formation Constants File", "Chemical Reaction Network File",
                   "Simulation Parameters File"]
        for s in l_paths:
            outfile.write(s+": "+kwargs[s]+"\n")
        outfile.write("CPU cores: "+str(kwargs['Cores'])+"\n")
        outfile.write("Starting Time: " + str(datetime.datetime.now()) +"\n")
        outfile.write("\n==========================\n")
        outfile.write("Chemical Reaction Network:\n")
        outfile.write("==========================\n")
        l_crn = ["Use Isomorphisms", "Reaction Energy Threshold (kcal/mol)",
                "Proton Difference Threshold", "Reference Reactions"]
        for s in l_crn:
            outfile.write(s+": "+str(kwargs[s])+"\n")

        outfile.write("\n==================\n")
        outfile.write("Speciation Models:\n")
        outfile.write("==================\n")
        l_spe = ["Ionic Strength (mol/L)", "Initial Concentration (mol/L)", "Range of pH", "Step of pH"
                 ,"Range of Simulated Models", "Formation Constants Referred to"]
        for s in l_spe:
            outfile.write(s+": "+str(kwargs[s])+"\n")
        outfile.write("\n==================\n")
        outfile.write("System Labels:\n")
        outfile.write("==================\n")
        outfile.write(str(kwargs["Labels"])+"\n")

def write_speciationparameters(kwargs):
    """
    Generates an output file with the variables and parameters chosen for a given simulation
    """


    with open(kwargs["Speciation Parameters File"], 'w') as outfile:


        outfile.write("=============================================================================================\n")
        outfile.write("·······································POMSimulator··········································\n")
        outfile.write("=============================================================================================\n")
        outfile.write("Information about the software's license, documentation and citation can be found at:   \n")
        outfile.write("https://gitlab.com/enricpp/pomsimulator_release                                         \n\n")
        outfile.write("Authors: Enric Petrus, Mireia Segado-Centellas, Carles Bo                                \n")
        outfile.write("Institute of Chemical Research of Catalonia (ICIQ) & Universitat Rovira i Virgili (URV) \n")
        outfile.write("Developer team: Enric Petrus, Jordi Buils, Diego Garay-Ruiz                                  \n")
        outfile.write("Contact Mail: enricpz@icloud.com                                                             \n")
        outfile.write("=============================================================================================\n\n")

        outfile.write("================\n")
        outfile.write("System Settings:\n")
        outfile.write("================\n")
        outfile.write("Formation constants file: %s \n" % kwargs["Formation Constants File"])
        outfile.write("Path to speciation output: %s \n" % kwargs["Path to Speciation Output"])
        outfile.write("Slope = %.4f Intercept = %.4f \n"%(kwargs["Scaling Slope"],kwargs["Scaling Intercept"]))
        outfile.write("Scaling type: %s \n"%kwargs["Scaling Type"])
        outfile.write("CPU cores: "+str(kwargs['Cores'])+"\n")
        outfile.write("Starting Time: " + str(datetime.datetime.now()) +"\n")
        outfile.write("\n==========================\n")


        l_spe = ["Initial Concentration (mol/L)", "Range of pH", "Step of pH","Number of Calculated Models", "Formation Constants Referred to"]
        for s in l_spe:
            outfile.write(s+": "+str(kwargs[s])+"\n")

def read_diagonal(path):
    """
    Reads existing isomorphic martrix and generates the corresponding numpt array
    Args:
        path: str, path to the subgraph isomorphism matrix
    Return:
        Diagonal: numpy array, isomorphic matrix
    """

    Diagonal = np.genfromtxt(path,delimiter=",")
    return Diagonal
