# Standard library imports
import networkx as nx
from itertools import product
import numpy as np
from multiprocessing import Pool
import matplotlib.pyplot as plt

# Local imports
from pomsimulator.modules.DataBase import Z_dict,Z_dict_inv, valence_dict

def Determine_IPA(atomic_number):
    """
    Determines which type of Isopolyoxoanions are present in the molecular set

    Args:
        atomic_number: list of integers, atomic numbers for a molecule <i>

    Returns:
        Z: integer, greatest atomic number for a molecule <i>
        Valence: integer, valence of the greatest atomic number for a molecule <i>
        Element: string, element
    """

    
    zlist = [obj for obj in atomic_number]
    Z = max(zlist)
    Element = Z_dict_inv[Z]
    Valence = valence_dict[Element]
    
    return Z, Valence, Element

def Molecule_Stoichiometry(G1_list, Z):
    """
    Computes the atom stoichiometry for a list of molecular graphs

    Args:
        G1_list: list of molecular graphs
        Z: integer, atomic number of the metal in the isopolyoxoanion

    Returns:
        stoich: list of lists, atom stoichiometries for every molecular graph
    """

    stoich = list()
    for i in range(len(G1_list)):
        if isinstance(Z, int):
            deleteme = list(nx.get_node_attributes(G1_list[i], 'Z').values())
            stoch_i = [deleteme.count(num) for num in [Z, 8, 1]]
            stoich.append(stoch_i)
        elif isinstance(Z, list):
            deleteme = list(nx.get_node_attributes(G1_list[i], 'Z').values())
            stoch_i = [deleteme.count(num) for num in Z+[8, 1]]
            stoich.append(stoch_i)

    return stoich

def Molecule_Charge(G1_list, Z, Zcharge):
    """
    Computes the molecular charges for a list of molecular graphs

    Args:
        G1_list: list of molecular graphs.
        Z: integer, atomic numbre of the metal in the isopolyoxoanion
        Zcharge: integer, charge of the metal in the isopolyoxoanion

    Returns:
        charge: list of lists, molecular charges for every molecular graph.
    """

    charge = list()
    for i in range(len(G1_list)):
        tmp = list(nx.get_node_attributes(G1_list[i], 'Z').values())
        stoichiometry = [tmp.count(num) for num in Z + [8, 1]]
        charge.append(sum([(v * valence) for v, valence in zip(stoichiometry, Zcharge + [-2, +1])]))
    return charge

def Molecule_to_Graph(idx, ploting=False, **kwargs): 
    """
    Converts a molecule to a molecular graph. The energy and name
    of the object are assigned as attributes.

    Args:
        idx: integer, index value for a molecule <i>.
        ploting: bolean, show a plot of the molecular graph
        kwards: dictionary, atomic numbers, bonds, label, and energies for a molecule <i>

    Returns:
        G: networkX object, molecular graph 

    """""

    atomic_number, connectivity, string = kwargs['Z'], kwargs['bonds'], kwargs['label']
    gibbs, enthalpy = kwargs['Gibbs'], kwargs['Enthalpy']
    
    Z, zmetal, metal = Determine_IPA(atomic_number)
    Z_list = [atomic_number.count(i) for i in [Z, 8, 1]]
    charge = sum([(v*valence) for v,valence in zip(Z_list, [zmetal, -2, +1])])
    G = nx.Graph(gibbs=gibbs, enthalpy=enthalpy, Label=idx, String=string)
    
    for i in range(len(atomic_number)):
        G.add_node(i + 1, Z=atomic_number[i])
    for bond in connectivity:
        G.add_edge(bond[0], bond[1])
    
    if ploting:
        nx.draw(G)
        plt.show()
    
    return G

def Molecule_to_Graph_from_molfile(idx, Z_list, bond_list, label):
    """
    Converts a molecule to a molecular graph. The energy and name
    of the object are assigned as attributes.

    Args:
        idx: Integer, index value for a molecule <i>.
        Z_list: Atomic numbers for a molecule
        bond_list: List of bonds in a molecules
        label: Name of molecule

    Returns:
        G: networkX object, molecular graph 

    """""

    G = nx.Graph(Index=idx, Label=label)

    for i,Z in enumerate(Z_list):
        G.add_node(i + 1, Z = Z)
    for bond in bond_list:
        G.add_edge(bond[0], bond[1])

    return G

def Create_Stoich(G1_labels):
    """
    Creates index groups according to the labeling of the output filess. This step is 
    necessary to ensure that the speciation models contain all the nuclearities.

    Side Note: needs to be refactored to create such objects from the stoich list
    
    Args:
        G1_labels: list of strings, names of the ADF output files
    
    Returns:
        unique_lab: list of strings, names of the ADF files for each nuclearity
        compounds_set: list of lists, number of molecules for each nuclearity

    """

    # Determine how many groups there are
    lab = tuple(s.split("-")[0] for s in G1_labels)
    unique_lab = sorted(list(set(lab)))   # CARE, set unsorts cause dictionary!
    lab_l = list(lab)
    cnt_lab = [lab_l.count(u) for u in unique_lab] # if u not in Monomers]
    compounds_set = [[] for _ in range(len(cnt_lab))]
    
    # Create a list with the stoichiometry groups
    acc = 0
    for ind, cnt in enumerate(cnt_lab):
        for i in range(cnt):
            num = acc + i
            compounds_set[ind].append(num)
        acc = num + 1
    return compounds_set, unique_lab

def Reaction_Type_IPA(G_list, M, ind, water, threshold, valid_cond, hydrated_species):
    """
    Reaction template for finding the chemical transformation of isopolyoxometalates. The
    user should tweak this function (or create a new one) depending on the needs of the system.
    It is meant as a modular template.

    Args:
        G_list: list of networkX objects, molecular graphs
        M: Metal name
        ind: integer, molecular graph index,
        water: dictionary, energies of water and its derivates
        threshold: integer, energy threshold to filter out reactions
        valid_cond: dictionary of conditions
        hydrated_species: list of hydrated species
    Returns:
        reac_e: reaction energy for a concrete reaction type
        reac_type: reaction type for a concrete reaction
    """
    
    ZM, ZO, ZH = [Z_dict[symb] for symb in [M, 'O', 'H']]
    if len(ind) == 2:  # Unimolecular Reactions
        p1, r1 = ind
        Z_p1 = list(nx.get_node_attributes(G_list[p1], 'Z').values())  # Atomic Numbers List
        Z_r1 = list(nx.get_node_attributes(G_list[r1], 'Z').values())
        v_p1 = [Z_p1.count(z) for z in [ZM, ZO, ZH]]  # Stoichiometry
        v_r1 = [Z_r1.count(z) for z in [ZM, ZO, ZH]]
        v_p2 = [(v_p1[m] - v_r1[m]) for m in range(3)]

        if v_p2 == [0, 0, 1]:  # M + H5O2 -> HM + 2H2O
            reac_e = G_list[p1].graph['gibbs'] + water['H4O2'] - (G_list[r1].graph['gibbs'] + water['H5O2'])
            return reac_e, 'P'

        elif v_p2 == [0, 1, 2]:  # MO4 + 2H2O -> MO4·2H2O
            reac_e = G_list[p1].graph['gibbs'] - (G_list[r1].graph['gibbs'] + 1 * water['H2O'])
            return reac_e, 'H2Ow1'
        elif v_p2 == [0, 2, 4]:  # MO4 + 2H2O -> MO4·2H2O
            reac_e = G_list[p1].graph['gibbs'] - (G_list[r1].graph['gibbs'] + 2 * water['H2O'])
            return reac_e, 'H2Ow2'
        elif v_p2 == [0, 1, 3]:  # MoO4 + H3O -> H3MoO5 (Acid hydration)
            reac_e = G_list[p1].graph['gibbs'] - (G_list[r1].graph['gibbs'] + water['H3O'])
            return reac_e, 'H3O'
        elif v_p2 == [0, 1, 1] and 'H6O3' in water.keys():  # M + 3H2O -> HM + H5O2 (Hydroxylation)
            reac_e = G_list[p1].graph['gibbs'] + water['H5O2'] - (G_list[r1].graph['gibbs'] + water['H6O3'])
            return reac_e, 'HO'
        else:
            return None

    elif len(ind) == 3:  # Bimolecular Reactions
        p1, r1, r2 = ind
        Z_p1 = list(nx.get_node_attributes(G_list[p1], 'Z').values())  # Atomic Numbers List
        Z_r1 = list(nx.get_node_attributes(G_list[r1], 'Z').values())
        Z_r2 = list(nx.get_node_attributes(G_list[r2], 'Z').values())
        v_p1 = [Z_p1.count(z) for z in [ZM, ZO, ZH]]  # Stoichiometry
        v_r1 = [Z_r1.count(z) for z in [ZM, ZO, ZH]]
        v_r2 = [Z_r2.count(z) for z in [ZM, ZO, ZH]]
        v_p2 = [(v_p1[m] - (v_r1[m] + v_r2[m])) for m in range(3)]

        # Handle hydrated species to check proton numbering if condition is selected
        if valid_cond.get("adjust_protons_hydration",False):
            vr_nw_list = []
            for vr in [v_r1,v_r2]:
                if tuple(vr) in hydrated_species.keys():
                    nwater = hydrated_species[tuple(vr)]
                    vr_nw = [vr[0],vr[1] - 1 * nwater,vr[2] - 2 * nwater]
                    vr_nw_list.append(vr_nw)
                else:
                    vr_nw_list.append(vr)
            v_r1_nw,v_r2_nw = vr_nw_list
            proton_dif = abs(v_r1_nw[2] - v_r2_nw[2])
        else:
            proton_dif = abs(v_r1[2] - v_r2[2])
            
        has_monomer = (v_r1[0] == 1 or v_r2[0] == 1)
        conditions_addition = True
        conditions_condensation = True
        force_condition = False

        if 'restrain_addition' in valid_cond.keys():
            val = valid_cond['restrain_addition']
            conditions_addition = conditions_addition and (v_r1[0] <= val or v_r2[0] <= val)

        if 'restrain_condensation' in valid_cond.keys():
            val = valid_cond['restrain_condensation']
            conditions_condensation = conditions_condensation and (v_r1[0] <= val or v_r2[0] <= val)

        if valid_cond.get('include_dimerization',False):
            conditions_addition = conditions_addition or (v_r1 == v_r2)
            conditions_condensation = conditions_condensation or (v_r1 == v_r2)

        if 'force_stoich' in valid_cond.keys():
            sel_stoich = valid_cond['force_stoich']
            force_condition = v_p1[0] in sel_stoich

        if 'proton_numb' in valid_cond.keys():
            conditions_addition = conditions_addition and (proton_dif <= valid_cond['proton_numb'])
            conditions_condensation = conditions_condensation and (proton_dif <= valid_cond['proton_numb'])

        # Classification
        if v_p2 == [0, 0, 0]: 
            """Addition Type: + HxMO4"""
            reac_e = G_list[p1].graph['gibbs'] - (G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_addition):
                return reac_e, 'A'

        elif v_p2 == [0, -1, -2]:
            """Condensation Type: MO4 + M2O7 --> M3O10 + H2O"""
            reac_e = G_list[p1].graph['gibbs'] + 1 * water['H2O'] - (
                        G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_condensation) or force_condition:
                return reac_e, 'Cw1'
        elif v_p2 == [0, -2, -4]:
            """Condensation Type: MO4 + M2O7 --> M3O9 + 2H2O"""
            reac_e = G_list[p1].graph['gibbs'] + 2 * water['H2O'] - (
                        G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_condensation) or force_condition:
                return reac_e, 'Cw2'
        elif v_p2 == [0, -3, -6]:
            """Condensation Type: M1O6 + M7O23 --> M8O26 + 3 H2O """
            reac_e = G_list[p1].graph['gibbs'] + 3 * water['H2O'] - (
                        G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_condensation) or force_condition:
                return reac_e, 'Cw3'
        elif v_p2 == [0, -4, -8]:
            """Condensation Type: M1O6 + M8O35 --> M9O37 + 4 H2O """
            reac_e = G_list[p1].graph['gibbs'] + 4 * water['H2O'] - (
                        G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_condensation) or force_condition:
                return reac_e, 'Cw4'
        elif v_p2 == [0, -10, -20]:
            reac_e = G_list[p1].graph['gibbs'] + 10 * water['H2O'] - (
                        G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_condensation) or force_condition:
                return reac_e, 'Cw10'
        else:
            return None

def Reaction_Type_HPA(G_list,XM, ind, water, threshold, valid_cond, hydrated_species):
    """
    Reaction template for finding the chemical transformation of heteropolyoxometalates. The
    user should tweak this function (or create a new one) depending on the needs of the system.
    It is meant as a modular template.

    Args:
        G_list: list of networkX objects, molecular graphs
        XM: Heteroatom and Metal name
        ind: integer, molecular graph index,
        water: dictionary, energies of water and its derivates
        threshold: integer, energy threshold to filter out reactions
        valid_cond: dictionary of conditions
        hydrated_species: list of hydrated species
    Returns:
        reac_e: reaction energy for a concrete reaction type
        reac_type: reaction type for a concrete reaction
    """

    X, M = XM.split('_')
    ZX, ZM, ZO, ZH = [Z_dict[symb] for symb in [X, M, 'O', 'H']]

    if len(ind) == 2:  # Unimolecular Reactions
        p1, r1 = ind
        Z_p1 = list(nx.get_node_attributes(G_list[p1], 'Z').values())  # Atomic Numbers List
        Z_r1 = list(nx.get_node_attributes(G_list[r1], 'Z').values())
        v_p1 = [Z_p1.count(z) for z in [ZX, ZM, ZO, ZH]]  # Stoichiometry
        v_r1 = [Z_r1.count(z) for z in [ZX, ZM, ZO, ZH]]
        v_p2 = [(v_p1[m] - v_r1[m]) for m in range(4)]

        if v_p2 == [0, 0, 0, 1]:  # M + H5O2 -> HM + 2H2O
            reac_e = G_list[p1].graph['gibbs'] + water['H4O2'] - (G_list[r1].graph['gibbs'] + water['H5O2'])
            return reac_e, 'P'
        elif v_p2 == [0, 0, 1, 2]:  # MO4 + 2H2O -> MO4·2H2O
            reac_e = G_list[p1].graph['gibbs'] - (G_list[r1].graph['gibbs'] + 1 * water['H2O'])
            return reac_e, 'H2Ow1'
        elif v_p2 == [0, 0, 2, 4]:  # MoO4 + 2H2O -> MoO4·2H2O
            reac_e = G_list[p1].graph['gibbs'] - (G_list[r1].graph['gibbs'] + 2 * water['H2O'])
            # if reac_e < max_perc:
            return reac_e, 'H2Ow2'
        elif v_p2 == [0, 0, 1, 3]:  # MoO4 + H3O -> H3MoO5 (Acid hydration)
            reac_e = G_list[p1].graph['gibbs'] - (G_list[r1].graph['gibbs'] + water['H3O'])
            return reac_e, 'H3O'
        elif v_p2 == [0, 0, 1, 1] and 'H6O3' in water.keys():  # M + 3H2O -> HM + H5O2 (Hydroxylation)
            reac_e = G_list[p1].graph['gibbs'] + water['H5O2'] - (G_list[r1].graph['gibbs'] + water['H6O3'])
            return reac_e, 'HO'
        else:
            return None

    elif len(ind) == 3:  # Bimolecular Reactions
        p1, r1, r2 = ind
        Z_p1 = list(nx.get_node_attributes(G_list[p1], 'Z').values())  # Atomic Numbers List
        Z_r1 = list(nx.get_node_attributes(G_list[r1], 'Z').values())
        Z_r2 = list(nx.get_node_attributes(G_list[r2], 'Z').values())
        v_p1 = [Z_p1.count(z) for z in [ZX, ZM, ZO, ZH]]  # Stoichiometry
        v_r1 = [Z_r1.count(z) for z in [ZX, ZM, ZO, ZH]]
        v_r2 = [Z_r2.count(z) for z in [ZX, ZM, ZO, ZH]]
        v_p2 = [(v_p1[m] - (v_r1[m] + v_r2[m])) for m in range(4)]


        # Handle hydrated species to check proton numbering if condition is selected
        if valid_cond.get("adjust_protons_hydration", False):
            vr_nw_list = []
            for vr in [v_r1, v_r2]:
                if tuple(vr) in hydrated_species.keys():
                    nwater = hydrated_species[tuple(vr)]
                    vr_nw = [vr[0], vr[1], vr[2] - 1 * nwater, vr[3] - 2 * nwater]
                    vr_nw_list.append(vr_nw)
                else:
                    vr_nw_list.append(vr)
            v_r1_nw, v_r2_nw = vr_nw_list
            proton_dif = abs(v_r1_nw[3] - v_r2_nw[3])
        else:
            proton_dif = abs(v_r1[3] - v_r2[3])

        has_metal_monomer = (v_r1[1] == 1 or v_r2[1] == 1)
        conditions_addition = True
        conditions_condensation = True
        force_condition = False

        if 'restrain_addition' in valid_cond.keys():
            val = valid_cond['restrain_addition']
            conditions_addition = conditions_addition and (v_r1[1] <= val or v_r2[1] <= val)

        if 'restrain_condensation' in valid_cond.keys():
            val = valid_cond['restrain_condensation']
            conditions_condensation = conditions_condensation and (v_r1[1] <= val or v_r2[1] <= val)

        if valid_cond.get('include_dimerization', False):
            conditions_addition = conditions_addition or (v_r1 == v_r2)
            conditions_condensation = conditions_condensation or (v_r1 == v_r2)

        if 'force_stoich' in valid_cond.keys():
            sel_stoich = valid_cond['force_stoich']
            force_condition = v_p1[1] in sel_stoich

        if 'proton_numb' in valid_cond.keys():
            conditions_addition = conditions_addition and (proton_dif <= valid_cond['proton_numb'])
            conditions_condensation = conditions_condensation and (proton_dif <= valid_cond['proton_numb'])

        # Classification
        if v_p2 == [0, 0, 0, 0]:
            """Addition Type: +HxMoO4"""
            reac_e = G_list[p1].graph['gibbs'] - (G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_addition):
                return reac_e, 'A'
        elif v_p2 == [0, 0, -1, -2]:
            """Condensation Type: +HxMoO4 + H2O"""
            reac_e = G_list[p1].graph['gibbs'] + 1 * water['H2O'] - (G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_condensation) or force_condition:
                return reac_e, 'Cw1'
        elif v_p2 == [0, 0, -2, -4]:
            """Condensation Type: +HxMoO4 + 2 H2O"""
            reac_e = G_list[p1].graph['gibbs'] + 2 * water['H2O'] - (G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_condensation) or force_condition:
                return reac_e, 'Cw2'
        elif v_p2 == [0, 0, -3, -6]:
            """Condensation Type: +HxMoO4 + 3 H2O"""
            reac_e = G_list[p1].graph['gibbs'] + 3 * water['H2O'] - (G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_condensation) or force_condition:
                return reac_e, 'Cw3'
        elif v_p2 == [0, 0, -4, -8]:
            """Condensation Type: +HxMoO4 + 4 H2O"""
            reac_e = G_list[p1].graph['gibbs'] + 4 * water['H2O'] - (G_list[r1].graph['gibbs'] + G_list[r2].graph['gibbs'])
            if (reac_e < threshold and conditions_condensation) or force_condition:
                return reac_e, 'Cw4'
        else:
            return None

def Isomorphism_to_ChemicalReactions(G1_list, np_IM, water, reference, POM, threshold, cond_dict_raw):
    """
    Returns a list of chemical reactions by processing a list of isomorphisms (i.e.,
    Isomorphic Matrix). Heuristics are applied for this particular case, where only
    acid-base, condensation and addition reactions exist. The transformation is
    carried out based on the stoichiometric difference between to isomorphic graphs.

    For example:

    i)   Graph_1 (H2MO4); Graph_2(H2M2O7)
    ii)  Graph_1 is subgraph of Graph_2, which means H2MO4 + X  -->  H2M2O7
    iii) Substracting the atoms: Graph_1 - Graph_2 = H0, M1, O3
    iv)  X = M1O3 which is: [MO4]2- and H2O

    Final result: H2MO4 + [MO4]2-  -->  H2M2O7 + H2O


    Args:
        G1_list: list of networkx objects, molecular graphs
        np_IM: list of list, isomorphic matrix
        water: dictionary, energies for water and its derivates
        reference: list of strings, reaction types considered in the network
        POM: string, type of isopolyoxometalate
        threshold: integer, energy threshold to filter out reactions
        cond_dict: dictionary of conditions

    Returns:
        Reac_idx: list of integers, combination of chemical reaction indexes.
        Reac_energy: list of floats, combination of chemical reaction energies.
        Reac_type: list of strings, combination of chemical reaction types.

    """

    if "_" in POM:
        func = Reaction_Type_HPA
    else:
        func = Reaction_Type_IPA
    conditions_list = ["proton_numb", "restrain_addition", "restrain_condensation", "include_dimerization",
                       "force_stoich", "adjust_protons_hydration"]

    # Format internal conditions
    valid_cond = {"proton_numb":int(cond_dict_raw["proton_numb"]),
                 "restrain_addition":int(cond_dict_raw["restrain_addition"]),
                 "restrain_condensation":int(cond_dict_raw["restrain_condensation"]),
                 "include_dimerization":cond_dict_raw.getboolean("include_dimerization"),
                 "adjust_protons_hydration":cond_dict_raw.getboolean("adjust_protons_hydration"),
                 "force_stoich":[int(item) if item else None for item in cond_dict_raw["force_stoich"].split(",")]
                 }
    for cond in cond_dict_raw.keys():
        if cond not in conditions_list:
            print("Condition %s is not supported. Valid conditions are %s" % (cond, ",".join(conditions_list)))

    reac_e_eq, found_reactions = list(), list()
    reac_info = list()
    N_Graphs = len(G1_list)
    print("IN PROCESS, PLEASE WAIT...")
    monomolec_comb = product(range(N_Graphs),repeat=2)
    bimolec_comb = product(range(N_Graphs),repeat=3)

    # Flag hydrated species so they are not considered when applying proton_numb filter
    Z_atoms = [Z_dict[symb] for symb in POM.split("_") + ['O', 'H']]
    hydrated_products = {}
    for i,j in monomolec_comb:
        obj = func(G1_list, POM, [i, j], water, threshold, valid_cond, hydrated_species=[])
        if obj == None or i == j:
            continue
        if np_IM[i][j] == 1:  # Unimolecular Reaction Found
            reac_g, reac_type = obj
            reac_info.append(([i, j], reac_g, reac_type))
            # get stoich to save hydrated species
            if reac_type in ["H2Ow1","H2Ow2"]:
                Z_hyd = list(nx.get_node_attributes(G1_list[i], 'Z').values())
                v_hyd = [Z_hyd.count(z) for z in Z_atoms]
                n_water = int(reac_type.replace("H2Ow",""))
                hydrated_products[tuple(v_hyd)] = n_water

    for i,j,k in bimolec_comb:
        obj = func(G1_list, POM, [i, j, k], water, threshold, valid_cond, hydrated_products)
        if obj == None:
            continue

        if [i, j, k] not in found_reactions:  # Bimolecular Reaction Found
            if np_IM[i][k] == 1 and np_IM[i][j] == 1:
                reac_g, reac_type = obj
                reac_info.append(([i, j, k], reac_g, reac_type))
                found_reactions.append([i, j, k])
                found_reactions.append([i, k, j])

    num = len(reference)
    Reac_idx, Reac_energy, Reac_type = [[] for _ in range(num)], [[] for _ in range(num)], [[] for _ in range(num)]
    ignored_reactions = list()
    for ind, e, t in reac_info:
        if t not in reference:
            ignored_reactions.append([ind,e,t])
            continue
        position = reference.index(t)
        Reac_idx[position].append(ind)
        Reac_energy[position].append(e)
        Reac_type[position].append(t)
    print("Number of Reactions", list(map(len, Reac_idx)))
    if ignored_reactions:
        print("%d reactions were not considered, as their types were not in reference" % len(ignored_reactions))
        print(ignored_reactions)
    return Reac_idx, Reac_energy, Reac_type

def _wrapper_isomorphism(G_pair):
    """
    Private function that Computes the isomorphism between a molecular graph (i) and a molecular graph (j). It 
    considers the atomic numbers as attributes when the isomorphism property is evaluated.

    Args:
        G_pair: list of two networkX objects, molecular graphs

    Returns:
        is_it_isomorphic: bolean, True is isomorphic
    """

    dict_bol = {True: 1, False: 0}
    if G_pair is None:
        return 0
    else:
        Gi, Gj = G_pair
        node_match = nx.algorithms.isomorphism.categorical_node_match('Z', 0)  # Matching Atomic Numbers
        GM = nx.algorithms.isomorphism.GraphMatcher(Gi, Gj, node_match=node_match)
        value = GM.subgraph_is_isomorphic()
    is_it_isomorphic = dict_bol[value]

    return is_it_isomorphic 

def _create_moleculargraph_withoutprotons(Gi):
    """
    Private function for removing protons of molecular graphs. It is employed to smooth the isomorphism search.
    This assumption is made based on the fact that protons are labile in aqueous solution.

    Args:
        Gi: networkX object, molecular graph with protons

    Returns:
        Gj: networkX object, molecular graph without protons
    """

    nodes, edges = Gi.nodes, Gi.edges
    z_dict = nx.get_node_attributes(Gi, 'Z')
    Gj = nx.Graph()
    Gj.graph = Gi.graph
    for key in z_dict:
        Z = z_dict[key]
        if Z != 1: # filtering protons
            Gj.add_node(key, Z=Z)
    for bond in edges:
        i, j = bond
        Z1, Z2 = z_dict[i], z_dict[j]
        if 1 not in [Z1, Z2]: #filtering bonds to protons
            Gj.add_edge(i, j)
    
    return Gj

def Molecular_Graphs_to_Isomorphic_Matrix(G1_list, diag, cores=1, verbose=True):
    """
    Computes the Isomorphic Matrix. To reduce the computational cost, we use the grid given by the diag matrix
    which already disregards the repeated reactions.

    Args:
        G1_list: list of networkX objects, molecular graphs
        diag: list of lists, initial guess of the isomorphic matrix
        cores: integer, number of cores used to run the isomorphic search
        verbose: bolean, shows an important warning

    Returns:
        np_IM: list of list, isomorphic matrix 

    """
    
    if verbose:
        print("WARNING: The bond connectivity from the QTAIM calculation might contain some artificial bonds.\n It is highly recommended to check the .mol files, and refine the connectivity if needed.")

    G_list = [_create_moleculargraph_withoutprotons(g) for g in G1_list]

    IM = list()
    args = zip([(gi, gj) if diag[i][j] == 1 else None for i,gi in enumerate(G_list) for j,gj in enumerate(G_list)])
    with Pool(cores) as ThreadPool:
        IM = IM + ThreadPool.starmap(_wrapper_isomorphism, args)
    IM_2d, tmp, j = [], [], 0
    num_molgraph = len(G_list)
    iterator = [diag[i][j] for i,gi in enumerate(G_list) for j,gj in enumerate(G_list)]
    for i,iso in enumerate(iterator):
        if j < num_molgraph:
            tmp.append(IM[i])
            j += 1
        if j == num_molgraph:
            IM_2d.append(tmp)
            j = 0
            tmp = list()

    np_IM = np.array(IM_2d)

    return np_IM

def sort_by_type(Reac_idx, Reac_energy, Reac_type, compounds_set, unique_labels, labels, exclude):
    """
    Sorts chemical reactions in a convenient manner so that speciation model can be build and solved in a sole
    <for> loop. This sorting eases the implementation later on of two chemical assumptions of POMSimulator:
    1) proton reactions are constant in each speciation model, 2) all nuclearities must be represented with one
    condensation/addition reaction

    Args:
        Reac_idx: list of integers, combination of chemical reaction indexes.
        Reac_energy: list of floats, combination of chemical reaction energies.
        Reac_type: list of strings, combination of chemical reaction types.
        unique_labels: list of strings, names of the ADF files for each nuclearity
        compounds_set: list of lists, number of molecules for each nuclearity
        labels: list of strings, name of the ADF outputs
        exclude: string, name of the reference compounds (used as reference in the formation constants)
    
    Returns:

        R_idx2: list of lists of integers, chemical reaction indexes organized by nuclearity
        R_ene2: list of lists of floats, chemical reaction energies organized by nuclearity
        R_type2: list of lists of strings, chemical reaction types organized by nuclearity
    """

    size = len(compounds_set)

    R_idx, R_ene, R_type = [[] for _ in range(size)], [[] for _ in range(size)], [[] for _ in range(size)]

    for ridx, eidx, tidx in zip(Reac_idx[1:], Reac_energy[1:], Reac_type[1:]):

        check = [[labels[obj2] for obj2 in obj] for obj in ridx]
        for ri, ei, ti in zip(ridx, eidx, tidx):
            prod = ri[0]

            for idx, c_set in enumerate(compounds_set):
                if prod in c_set:
                    R_idx[idx].append(ri)
                    R_ene[idx].append(ei)
                    R_type[idx].append(ti)

    compounds_set = tuple((tuple(obj) for obj in compounds_set))

    R_idx2, R_ene2, R_type2 = list(), list(), list()

    # acc = 0
    for a,b,c,d in zip(unique_labels, R_idx, R_ene, R_type):
        if a in exclude or len(b) == 0:
            print(a," was excluded")
        else:
            R_idx2.append(b)
            R_ene2.append(c)
            R_type2.append(d)
        # acc += 1
    print("Lenght unique labels", len(unique_labels))
    print("Number of Reactions for each Nuclearity",  dict([(a,b) for a,b in zip(unique_labels, list(map(len, R_idx)))])) #len(R_idx2), len(R_ene2), len(R_type2),
    return R_idx2, R_ene2, R_type2


