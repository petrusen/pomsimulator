"""
   Database file with the experimental data employed in the calibration step.
   Stoichiometric formulas are also depicted in this file.

"""

# EXPERIMENTAL FORMATION CONSTANTS for Tungsten polyoxometalates reported by Rozantzev and Sazonova, Russ. J. Cord. Chem. 2005, 31 ,552
Rosantsev_W12_I01_05I = {"W06O22-2H": 53.68, 
                         "W07O24-1H": 76.59, 
                         "W12O40-2H": 149.59,
                         "W10O32-0H": 129.63}
Pettersson_3I = {"P00Mo01O04-0H":0,
                 "P01Mo00O04-0H":0,
                 "P01Mo09O34-2H":104.9,
                 "P01Mo09O34-1H":102.0,
                 "P01Mo11O39-0H":118.7,
                 "P01Mo12O40-0H":139.7,
                 "P02Mo05O23-0H": 61.97,
                 "P02Mo05O23-1H": 67.07,
                 "P02Mo05O23-2H": 70.86}

experimental_constants = {"W12_Rosantsev_I01_05":Rosantsev_W12_I01_05I,
    "PMo12_Petterson_I3":Pettersson_3I}

# COLORS FOR THE PLOTS

Col_Dict_PMo = {"P00Mo01O04-0H": "#393b79ff",
            "P00Mo01O04-1H": "#393b79ff",
            "P00Mo01O04-2H": "#393b79ff",
            "P00Mo01O04-3H": "#393b79ff",
            "P00Mo01O06-4H": "#5254a3ff",
            "P00Mo01O06-5H": "#5254a3ff",
            "P00Mo01O06-6H": "#5254a3ff",
            "P00Mo01O06-7H": "#5254a3ff",
            "P00Mo01O06-8H": "#5254a3ff",
            "P00Mo02O06-0H": "#6b6ecfff",
            "P00Mo02O06-1H": "#6b6ecfff",
            "P00Mo02O06-2H": "#6b6ecfff",
            "P00Mo02O07-0H": "#9294daff",
            "P00Mo02O07-1H": "#9294daff",
            "P00Mo02O07-2H": "#9294daff",
            "P00Mo02O08-0H": "#c3c3e3ff",
            "P00Mo02O08-1H": "#c3c3e3ff",
            "P00Mo02O08-2H": "#c3c3e3ff",
            "P00Mo03O09-0H": "#637939ff",
            "P00Mo03O09-1H": "#637939ff",
            "P00Mo03O09-2H": "#637939ff",
            "P00Mo03O10-0H": "#8ca252ff",
            "P00Mo03O10-1H": "#8ca252ff",
            "P00Mo03O10-2H": "#8ca252ff",
            "P00Mo03O11-0H": "#b5cf6bff",
            "P00Mo03O11-1H": "#b5cf6bff",
            "P00Mo03O11-2H": "#b5cf6bff",
            "P00Mo04O13-0H": "#cedb9cff",
            "P00Mo04O13-1H": "#cedb9cff",
            "P00Mo04O13-2H": "#cedb9cff",
            "P00Mo05O16-0H": "#8c6d31ff",
            "P00Mo05O16-1H": "#8c6d31ff",
            "P00Mo05O16-2H": "#8c6d31ff",
            "P00Mo05O17-0H": "#bd9e39ff",
            "P00Mo05O17-1H": "#bd9e39ff",
            "P00Mo05O17-2H": "#bd9e39ff",
            "P00Mo06O20-0H": "#e7ba52ff",
            "P00Mo06O20-1H": "#e7ba52ff",
            "P00Mo06O21-0H": "#e7cb94ff",
            "P00Mo06O21-1H": "#e7cb94ff",
            "P00Mo06O21-2H": "#e7cb94ff",
            "P01Mo00O04-0H": "#843c39ff",
            "P01Mo00O04-1H": "#ad494aff",
            "P01Mo00O04-2H": "#d6616bff",
            "P01Mo00O04-3H": "#e7969cff",
            "P01Mo03O13-0H": "#7b4173ff",
            "P01Mo03O13-1H": "#a55194ff",
            "P01Mo03O13-2H": "#ce6dbdff",
            "P01Mo03O13-3H": "#de9ed6ff",
            "P01Mo05O19-0H": "#82f7e3ff",
            "P01Mo05O19-1H": "#82f7e3ff",
            "P01Mo05O19-2H": "#b9eae2ff",
            "P01Mo06O22-0H": "#555555ff",
            "P01Mo06O22-1H": "#808080ff",
            "P01Mo06O22-2H": "#aeaeaeff",
            "P01Mo06O22-3H": "#e2e2e2ff",
            "P01Mo09O31-0H": "#e48308ff",
            "P01Mo09O31-1H": "#e48308ff",
            "P01Mo09O31-2H": "#f8a43dff",
            "P01Mo09O31-3H": "#f8a43dff",
            "P01Mo09O34-0H": "#fbca8bff",
            "P01Mo09O34-1H": "#fbca8bff",
            "P01Mo09O34-2H": "#fde3c1ff",
            "P01Mo09O34-3H": "#fde3c1ff",
            "P01Mo09O34-5H": "#fde3c1ff",
            "P01Mo09O34-6H": "#fde3c1ff",
            "P01Mo11O39-0H": "#40af47ff",
            "P01Mo11O39-1H": "#5acf60ff",
            "P01Mo11O39-2H": "#7fe884ff",
            "P01Mo11O39-3H": "#b3f1b6ff",
            "P01Mo11O39-4H": "#b3f1b6ff",
            "P01Mo12O40-0H": "#ae266bff",
            "P01Mo12O40-1H": "#d84c93ff",
            "P01Mo12O40-2H": "#e37db1ff",
            "P02Mo05O23-0H": "#04b899ff",
            "P02Mo05O23-1H": "#33dbc0ff",
            "P02Mo05O23-2H": "#33dbc0ff",
            'P00Mo00O00-0H': "#000002ff"}

color_dictionaries = {"PMo_Keggin": Col_Dict_PMo}

allowed_scaling_modes = ["best_rmse","average","medians","universal"]

# Periodic table mapping
Z_dict = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12,
          'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23,
          'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
          'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45,
          'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56,
          'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
          'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
          'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89,
          'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
          'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds ': 110,
          'Rg ': 111, 'Cn ': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}

# Periodic table reverse mapping
Z_dict_inv = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 
              13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V',
              24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se',
              35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh',
              46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba',
              57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 
              68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
              79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac',
              90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
              101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds ',
              111: 'Rg ', 112: 'Cn ', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'}

# Oxidation states for common atoms present in polyoxometalates
valence_dict = {'W': 6, 'O': -2, 'H':1, 'Mo':6, 'P':5}

# Dictionary for mapping reaction type to its proper written from
stringreac_dict = {'P' : 'Acid Base: {R1} + H+ --> {P}    G={G}\n',
               'HO': 'Hydrox.: {R1} + H6O3 --> {P} + H5O2+  G={G}\n',
               'H2Ow1': 'Hydr. : {R1} + 1 H2O --> {P}    G={G}\n',
               'H2Ow2': 'Hydr. : {R1} + 2 H2O --> {P}    G={G}\n',
               'H3O': 'Acid Hydr. : {R1} + H3O+ --> {P}    G={G}\n',
               'Cw1': 'Cond. : {R1} + {R2} --> {P} + 1 H2O   G={G}\n',
               'Cw2': 'Cond. : {R1} + {R2} --> {P} + 2 H2O   G={G}\n',
               'Cw3': 'Cond. : {R1} + {R2} --> {P} + 3 H2O   G={G}\n',
               'Cw4': 'Cond. : {R1} + {R2} --> {P} + 4 H2O   G={G}\n',
               'Cw10': 'Cond. : {R1} + {R2} --> {P} + 10 H2O   G={G}\n',
               'A': 'Addition: {R1} + {R2} --> {P}  G={G}\n'}

# Dictionary for mapping reaction type to its equilibrium equation
equation_dict = {
    "P": "p[{a}] * h2o - exp(-{b} / (R * T)) * p[{c}] * c_H",
    "H2Ow1": "p[{a}] - exp(-{b} / (R * T)) * p[{c}] * h2o",
    "H2Ow2": "p[{a}] - exp(-{b} / (R * T)) * p[{c}] * h2o ** 2",
    "HO": "p[{a}] * c_H - exp(-{b} / (R * T)) * p[{c}] * h2o",
    "Cw1": "p[{a}] * h2o ** 1 - exp(-{b} / (R * T)) * p[{c}] * p[{d}]",
    "Cw2": "p[{a}] * h2o ** 2 - exp(-{b} / (R * T)) * p[{c}] * p[{d}]",
    "Cw3": "p[{a}] * h2o ** 3 - exp(-{b} / (R * T)) * p[{c}] * p[{d}]",
    "Cw4": "p[{a}] * h2o ** 4 - exp(-{b} / (R * T)) * p[{c}] * p[{d}]",
    "Cw10": "p[{a}] * h2o ** 10 - exp(-{b} / (R * T)) * p[{c}] * p[{d}]",
    "A": "p[{a}] - exp(-{b} / (R * T)) * p[{c}] * p[{d}]",
    "H3O": "p[{a}] - exp(-{b} / (R * T)) * p[{c}] * c_H" }

reaction_references = ["P", "H2Ow1", "H2Ow2", "Cw1", "Cw2", "Cw3", "Cw4", "A", "HO", "H3O"]

clustering_features = ["width","height","pos","area"]

def get_stoichiometry_K(ntotal, ener, idx, type, verbose=False):
    """
    Wrapper function to construct objects: one 2D-matrix containing all the reactions in the speciation
    model, and second, a simple array of the equilibrium constants.
    """
    R, T = 8.314 * 0.001 * (1 / 4.18), 298.15
    if type == "P":  # (ntotal-1= h2o) and (ntotal = h3o+=
        p1, r1 = idx
        tdic = {r1: -1, ntotal: -1, p1: 1, ntotal-1: 1}
        row = [tdic[i] if i in tdic.keys() else 0 for i in range(1, ntotal+1)], np.exp(-ener/(R*T))
    elif type == "Cw1":  # (ntotal-1= h2o) and (ntotal = h3o+=
        p1, r1, r2 = idx
        tdic = {r1: -1, r2: -1, ntotal-1: 1, p1: 1} if r1 != r2 else {r1: -2, ntotal-1: 1, p1: 1}
        row = [tdic[i] if i in tdic.keys() else 0 for i in range(1, ntotal+1)], np.exp(-ener/(R*T))
    elif type == "Cw2":  # (ntotal-1= h2o) and (ntotal = h3o+=
        p1, r1, r2 = idx
        tdic = {r1: -1, r2: -1, ntotal-1: 2, p1: 1} if r1 != r2 else {r1: -2, ntotal-1: 2, p1: 1}
        row = [tdic[i] if i in tdic.keys() else 0 for i in range(1, ntotal+1)], np.exp(-ener/(R*T))
    elif type == "Cw3":  # (ntotal-1= h2o) and (ntotal = h3o+=
        p1, r1, r2 = idx
        tdic = {r1: -1, r2: -1, ntotal-1: 3, p1: 1} if r1 != r2 else {r1: -2, ntotal-1: 3, p1: 1}
        row = [tdic[i] if i in tdic.keys() else 0 for i in range(1, ntotal+1)], np.exp(-ener/(R*T))
    elif type == "Cw4":  # (ntotal-1= h2o) and (ntotal = h3o+=
        p1, r1, r2 = idx
        tdic = {r1: -1, r2: -1, ntotal-1: 4, p1: 1} if r1 != r2 else {r1: -2, ntotal-1: 4, p1: 1}
        row = [tdic[i] if i in tdic.keys() else 0 for i in range(1, ntotal+1)], np.exp(-ener/(R*T))
    elif type == "A":  # (ntotal-1= h2o) and (ntotal = h3o+=
        p1, r1, r2 = idx
        tdic = {r1: -1, r2: -1, p1: 1} if r1 != r2 else {r1: -2, p1: 1}
        row = [tdic[i] if i in tdic.keys() else 0 for i in range(1, ntotal+1)], np.exp(-ener/(R*T))
    elif type == "H2Ow1":  # (ntotal-1= h2o) and (ntotal = h3o+=
        p1, r1 = idx
        tdic = {r1: -1, ntotal-1: -1, p1: 1}
        row = [tdic[i] if i in tdic.keys() else 0 for i in range(1, ntotal+1)], np.exp(-ener/(R*T))
    elif type == "H2Ow2":  # (ntotal-1= h2o) and (ntotal = h3o+=
        p1, r1 = idx
        tdic = {r1: -1, ntotal-1: -2, p1: 1}
        row = [tdic[i] if i in tdic.keys() else 0 for i in range(1, ntotal+1)], np.exp(-ener/(R*T))
    else:
        print(type, "=============================================================================")

    if verbose: print("row", type, row[0])
    return row

# Dictionary for mapping reaction types to their molecularity
molecularity_dict = {'Cw10': 2, 'Cw2': 2, 'Cw3': 2, 'Cw1': 2, 'Cw4': 2, 'A': 2, 'P' : 1, 'HO': 1, 'H2Ow1': 1, 'H2Ow2': 1, 'H3O' : 1}

# List of strings needed to write outfiles for simulations and speciations
simulation_parameters_strings = ['ADF Folder', 'MOL Folder', "Formation Constants File", "Chemical Reaction Network File",
            "Simulation Parameters File", "Cores", "Use Isomorphisms", "Reaction Energy Threshold (kcal/mol)",
            "Proton Difference Threshold", "Reference Reactions", "Ionic Strength (mol/L)",
            "Initial Concentration (mol/L)",  "Range of pH", "Step of pH", "Range of Simulated Models",
             "Formation Constants Referred to","Labels","Start date", "End date", "Execution time"]
speciation_parameters_strings = [ "Speciation Parameters File","Formation Constants File","Scaling Slope","Scaling Intercept","Scaling Type", "Cores",
            "Initial Concentration (mol/L)",  "Range of pH", "Step of pH", "Number of Calculated Models","Labels", "Formation Constants Referred to",
                                  "Path to Speciation Output", "Start date", "End date", "Execution time"]

# Universal scaling methodology
universal_slope = 0.29
"""                   Q3    mean   range  Indep."""
mlr_coefficients = [0.195, -0.216, 0.070, 12.20]

element_masses = {
    'H': 1.0079,
    'He': 4.0026,
    'Li': 6.941,
    'Be': 9.0122,
    'B': 10.81,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    'F': 18.998,
    'Ne': 20.180,
    'Na': 22.990,
    'Mg': 24.305,
    'Al': 26.982,
    'Si': 28.085,
    'P': 30.974,
    'S': 32.06,
    'Cl': 35.45,
    'Ar': 39.948,
    'K': 39.098,
    'Ca': 40.078,
    'Sc': 44.956,
    'Ti': 47.867,
    'V': 50.942,
    'Cr': 51.996,
    'Mn': 54.938,
    'Fe': 55.845,
    'Co': 58.933,
    'Ni': 58.693,
    'Cu': 63.546,
    'Zn': 65.38,
    'Ga': 69.723,
    'Ge': 72.63,
    'As': 74.922,
    'Se': 78.971,
    'Br': 79.904,
    'Kr': 83.798,
    'Rb': 85.468,
    'Sr': 87.62,
    'Y': 88.906,
    'Zr': 91.224,
    'Nb': 92.906,
    'Mo': 95.95,
    'Tc': 98.0,
    'Ru': 101.07,
    'Rh': 102.91,
    'Pd': 106.42,
    'Ag': 107.87,
    'Cd': 112.41,
    'In': 114.82,
    'Sn': 118.71,
    'Sb': 121.76,
    'Te': 127.60,
    'I': 126.90,
    'Xe': 131.29,
    'Cs': 132.91,
    'Ba': 137.33,
    'La': 138.91,
    'Ce': 140.12,
    'Pr': 140.91,
    'Nd': 144.24,
    'Pm': 145.0,
    'Sm': 150.36,
    'Eu': 151.96,
    'Gd': 157.25,
    'Tb': 158.93,
    'Dy': 162.50,
    'Ho': 164.93,
    'Er': 167.26,
    'Tm': 168.93,
    'Yb': 173.05,
    'Lu': 174.97,
    'Hf': 178.49,
    'Ta': 180.95,
    'W': 183.84,
    'Re': 186.21,
    'Os': 190.23,
    'Ir': 192.22,
    'Pt': 195.08,
    'Au': 196.97,
    'Hg': 200.59,
    'Tl': 204.38,
    'Pb': 207.2,
    'Bi': 208.98,
    'Po': 209.0,
    'At': 210.0,
    'Rn': 222.0,
    'Fr': 223.0,
    'Ra': 226.0,
    'Ac': 227.0,
    'Th': 232.04,
    'Pa': 231.04,
    'U': 238.03,
    'Np': 237.0,
    'Pu': 244.0,
    'Am': 243.0,
    'Cm': 247.0,
    'Bk': 247.0,
    'Cf': 251.0,
    'Es': 252.0,
    'Fm': 257.0,
    'Md': 258.0,
    'No': 259.0,
    'Lr': 262.0,
    'Rf': 267.0,
    'Db': 270.0,
    'Sg': 271.0,
    'Bh': 270.0,
    'Hs': 277.0,
    'Mt': 276.0,
    'Ds': 281.0,
    'Rg': 280.0,
    'Cn': 285.0,
    'Nh': 284.0,
    'Fl': 289.0,
    'Mc': 288.0,
    'Lv': 293.0,
    'Ts': 294.0,
    'Og': 294.0
}