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

# Labels extracted from molecular set, following the guidelines in the manual for species labeling.
Labels_W  = ['W01O04-0H', 'W01O04-1H', 'W01O04-2H', 
            'W01O06-6H', 'W01O06-7H', 'W01O06-8H',
            'W02O07-0H', 'W02O07-1H', 'W02O07-2H',
            'W03O10-0H', 'W03O10-1H', 'W03O10-2H',
            'W03O11-0H', 'W03O11-1H', 'W03O11-2H',
            'W04O13-0H', 'W04O13-1H', 'W04O13-2H',
            'W04O15-0H', 'W04O15-1H', 'W04O15-2H',
            'W05O16-0H', 'W05O16-1H', 'W05O16-2H',
            'W05O17-0H', 'W05O17-1H', 'W05O17-2H',
            'W05O19-0H', 'W05O19-1H', 'W05O19-2H',
            'W06O20-0H', 'W06O20-1H', 'W06O20-2H',
            'W06O22-0H', 'W06O22-1H', 'W06O22-2H',
            'W07O24-0H', 'W07O24-1H', 'W07O24-2H',
            'W10O32-0H', 'W10O32-1H', 'W10O32-2H',
            'W12O40-0H', 'W12O40-1H', 'W12O40-2H']


Labels_PMo = ['P00Mo01O04-0H', 'P00Mo01O04-1H', 'P00Mo01O04-2H',
              'P00Mo02O07-0H', 'P00Mo02O07-1H', 'P00Mo02O07-2H',
              'P00Mo02O08-1H', 'P00Mo02O08-2H',
              'P00Mo03O09-0H', 'P00Mo03O09-1H', 'P00Mo03O09-2H',
              'P00Mo03O10-0H', 'P00Mo03O10-1H', 'P00Mo03O10-2H',
              'P00Mo03O11-1H', 'P00Mo03O11-2H',
              'P00Mo06O21-1H', 'P00Mo06O21-2H',
              'P01Mo00O04-0H', 'P01Mo00O04-1H', 'P01Mo00O04-2H', 'P01Mo00O04-3H',
              'P01Mo03O13-0H', 'P01Mo03O13-1H', 'P01Mo03O13-2H', 'P01Mo03O13-3H',
              'P01Mo05O19-0H', 'P01Mo05O19-1H', 'P01Mo05O19-2H',
              'P01Mo06O22-0H', 'P01Mo06O22-1H',
              'P01Mo09O31-0H', 'P01Mo09O31-1H', 'P01Mo09O31-2H', 'P01Mo09O31-3H',
              'P01Mo09O34-0H', 'P01Mo09O34-1H', 'P01Mo09O34-2H', 'P01Mo09O34-3H',
              'P01Mo11O39-0H', 'P01Mo11O39-1H', 'P01Mo11O39-2H', 'P01Mo11O39-3H',
              'P01Mo12O40-0H', 'P01Mo12O40-1H', 'P01Mo12O40-2H',
              'P02Mo05O23-0H', 'P02Mo05O23-1H', 'P02Mo05O23-2H']

# Labels corresponding to species present in experimental studies
Labels_W_good = ['W01O04-0H','W06O22-2H', 'W07O24-1H', 'W12O40-2H', 'W10O32-0H']

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

# Dictionary for mapping reaction types to their molecularity
molecularity_dict = {'Cw10': 2, 'Cw2': 2, 'Cw3': 2, 'Cw1': 2, 'Cw4': 2, 'A': 2, 'P' : 1, 'HO': 1, 'H2Ow1': 1, 'H2Ow2': 1, 'H3O' : 1}

# List of strings needed to write outfiles for simulations and speciations
simulation_parameters_strings = ['ADF Folder', 'MOL Folder', "Formation Constants File", "Chemical Reaction Network File",
            "Simulation Parameters File", "Cores", "Use Isomorphisms", "Reaction Energy Threshold (kcal/mol)",
            "Proton Difference Threshold", "Reference Reactions", "Ionic Strength (mol/L)",
            "Initial Concentration (mol/L)",  "Range of pH", "Step of pH", "Range of Simulated Models", "Formation Constants Referred to","Labels"]
speciation_parameters_strings = [ "Speciation Parameters File","Formation Constants File","Scaling Slope","Scaling Intercept","Scaling Type", "Cores",
            "Initial Concentration (mol/L)",  "Range of pH", "Step of pH", "Number of Calculated Models","Labels", "Formation Constants Referred to",
                                  "Path to Speciation Output"]

# Universal scaling methodology
universal_slope = 0.29
"""                   Q3    mean   range  Indep."""
mlr_coefficients = [0.195, -0.216, 0.070, 12.20]

# Coloring
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
            'P00Mo00O00-0H': "#000000ff"}
