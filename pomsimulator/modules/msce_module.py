# Standard library imports
import numpy as np
from math import sqrt, log10, exp
from scipy.stats import linregress
from scipy.optimize import root,leastsq
from random import uniform
from sklearn.metrics import mean_squared_error
from itertools import repeat
import warnings
import os

# Local imports
from pomsimulator.modules.text_module import Lab_to_stoich
from pomsimulator.modules.DataBase import equation_dict,molecularity_dict


## Parameters to handle proper multithreading with numpy
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

warnings.filterwarnings("ignore")


## Management of parallelization
def starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    """Creates a wrapper around pool.starmap that also accepts an iterator
    of over kwargs dictionaries."""

    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(apply_args_and_kwargs, args_for_starmap)

def apply_args_and_kwargs(fn, args, kwargs):
    return fn(*args, **kwargs)

## Calculation of rate constants
def Speciation_from_Equilibrium(idx_var, e_var, type_var, idx_ctt=None, e_ctt=None, type_ctt=None, z_ctt=None, v_ctt=None,
                                ref_idx=None, pH_grid=None, init_guess=None, I=None, C=0.005, temp=298.15, solver='hybr', threshold=None):
    """
    Sets multi-species chemical equilibrium provided that the reactions are
    acid-base, condensations or addition reactions. The system of non-lineal
    equations are minimized using Powell's algorithm[1].

    The parameters have been fixed to the following default values:

        ·precision=5; the number of pH values at which MSCE are solved.
        ·R=0.00198; atmospheric pressure.
        ·[H2O]=1M; concentration of water.

    [1] Powell, M. J. D. An Efficient Method for Finding the Minimum of a Function of Several
    Variables without Calculating Derivatives. The Computer Journal. 1964, p 155.

    Args:
        idx_var: list of integers, combination of chemical reaction indexes.
        e_var: list of floats, combination of chemical reaction energies.
        type_var: list of strings, combination of chemical reaction types.
        idx_ctt, e_ctt, type_ctt: same as three first args, but for the constant reactions (acid_base).
        z_ctt: list of integers, charges of all compounds in the dataset.
        v_ctt: list of lists of integers, stoichiometric indices for chemical compounds in the dataset (metal, oxygen, hydrogen).
        ref_idx: integer, index of the reference species in compound-related lists.
        pH_grid: array of floats, all pH values to solve equations for.
        init_guess: array of floats, initial guess for the equation solver.
        I: float, ionic strength for the simulation.
        C: float, total metal concentration in the simulation.
        temp: float, temperature in Kelvin
        solver: string, method to be used for nonlinear system resolution from SciPy.
        threshold: float, controls maximum accepted rmse for a solution.
        verbose: boolean controlling the verbosity of the output.

    Returns:
        Kf_dft: list of floats, computed formation constants for all compounds in the reaction network.

    """

    solved_pH_val, solved_activity_val, = list(), list()
    reac_e_eq = e_ctt + list(e_var)
    reac_idx = idx_ctt + list(idx_var)
    reac_type = type_ctt + list(type_var)
    M_ratio = [v[0] for v in v_ctt]

    rmse_l = list()

    ### The system of equations is defined through strings, where values for indices and energies are substituted
    ### These strings are then precompiled as Python objects to improve performance and evaluated in the solving function
    compiled_eqs_list = list()
    for ener, idx, type in zip(reac_e_eq, reac_idx, reac_type):
        if molecularity_dict[type] ==1:
            eq = equation_dict[type].format(a=idx[0],b=ener,c=idx[1])
        else:
            eq = equation_dict[type].format(a=idx[0],b=ener,c=idx[1],d=idx[2])
        compiled_eqs_list.append(compile(eq,'<string>', 'eval'))

    mass_eq = ''.join([str(M_ratio[i]) + " * p[" + str(i) + "] + " for i in range(len(reac_e_eq) + 1)]) + str(-C)
    compiled_eqs_list.append(compile(mass_eq, '<string>', 'eval'))

    for idx, pH in enumerate(pH_grid):
        def non_lineal_sys(p):
            """
            Encapsulated function which iteratively sets and solves different systems
            of non-linear equations.

            Args:
                p: vector of concentrations, solved iteratively

            Returns:
                sys_eq. list of floats, solution of the equations depicted in the system.
            """
            R, T, c_H, h2o = 8.314 * 0.001 * (1 / 4.18), temp, 10 ** (- pH), 1
            sys_eq =  list()

            for k,eq in enumerate(compiled_eqs_list):  # sys_eq list comprehension does not work
                sys_eq.append(eval(eq))
            sys_eq = np.array(sys_eq)
            return sys_eq

        conc_i = root(non_lineal_sys, init_guess, method=solver)

        if conc_i.success == True:
            activities = Activity(conc_i.x, z_ctt, I)
            if np.isnan(activities).any() or np.isinf(activities).any():
                pass
            else:
                r2_i, rmse, mae, y_ls_values = Least_Squared(non_lineal_sys, init_guess, activities)
                rmse_acc = 0
                if rmse < (C * threshold):
                    solved_activity_val.append(activities)
                    solved_pH_val.append(pH)
                    rmse_acc = rmse_acc + rmse
                rmse_l.append(rmse_acc)

    if len(solved_activity_val) > 0:
        solved_activity_val_T = np.array(solved_activity_val).T
    else:
        solved_activity_val_T = list()

    Kf_dft = screen_log_Kf(solved_activity_val_T, solved_pH_val, v_ctt, ref_idx=ref_idx)

    return Kf_dft

def Speciation_from_Equilibrium_bimetal(idx_var, e_var, type_var, idx_ctt=None, e_ctt=None, type_ctt=None, z_ctt=None, v_ctt=None,
                                ref_idx=None, pH_grid=None, init_guess=None, I=None,
                                C_X=0.005, C_M=0.005, temp=298.15, solver='hybr', threshold=None):
    """
    Sets multi-species chemical equilibrium provided that the reactions are
    acid-base, condensations or addition reactions. The system of non-lineal
    equations are minimized using Powell's algorithm[1].

    The parameters have been fixed to the following default values:

        ·precision=5; the number of pH values at which MSCE are solved.
        ·R=0.00198; atmospheric pressure.
        ·[H2O]=1M; concentration of water.

    [1] Powell, M. J. D. An Efficient Method for Finding the Minimum of a Function of Several
    Variables without Calculating Derivatives. The Computer Journal. 1964, p 155.

    Args:
        idx_var: list of integers, combination of chemical reaction indexes.
        e_var: list of floats, combination of chemical reaction energies.
        type_var: list of strings, combination of chemical reaction types.
        idx_ctt, e_ctt, type_ctt: same as three first args, but for the constant reactions (acid_base).
        z_ctt: list of integers, charges of all compounds in the dataset.
        v_ctt: list of lists of integers, stoichiometric indices for chemical compounds in the dataset (metal, oxygen, hydrogen).
        ref_idx: integer list, index of the reference species in compound-related lists.
        pH_grid: array of floats, all pH values to solve equations for.
        init_guess: array of floats, initial guess for the equation solver.
        I: float, ionic strength for the simulation.
        C_X, C_M: float, total heteroatom and metal concentrations in the simulation.
        temp: float, temperature in Kelvin
        solver: string, method to be used for nonlinear system resolution from SciPy.
        threshold: float, controls maximum accepted rmse for a solution.
        verbose: boolean controlling the verbosity of the output.

    Returns:
        Kf_dft: list of floats, computed formation constants for all compounds in the reaction network.

        """

    solved_pH_val, solved_activity_val, rmse_l = list(), list(), list()
    reac_e_eq = e_ctt + list(e_var)
    reac_idx = idx_ctt + list(idx_var)
    reac_type = type_ctt + list(type_var)

    X_Ratio = [v[0] for v in v_ctt]
    M_Ratio = [v[1] for v in v_ctt]

    equations_list = list()
    compiled_eqs_list = list()
    for ener, idx, type in zip(reac_e_eq, reac_idx, reac_type):
        if type in ["P","H2Ow1","H2Ow2","H3O","HO"]:
            eq = equation_dict[type].format(a=idx[0],b=ener,c=idx[1])
        else:
            eq = equation_dict[type].format(a=idx[0],b=ener,c=idx[1],d=idx[2])

        equations_list.append(eq)
        compiled_eqs_list.append(compile(eq,'<string>', 'eval'))

    mass_eq1 = ''.join([str(X_Ratio[i]) + " * p[" + str(i) + "] + " for i in range(len(X_Ratio))]) + str(-C_X)
    mass_eq2 = ''.join([str(M_Ratio[i]) + " * p[" + str(i) + "] + " for i in range(len(M_Ratio))]) + str(-C_M)

    equations_list.append(mass_eq1)
    equations_list.append(mass_eq2)
    compiled_eqs_list.append(compile(mass_eq1,'<string>','eval'))
    compiled_eqs_list.append(compile(mass_eq2, '<string>', 'eval'))

    for idx, pH in enumerate(pH_grid):
        def non_lineal_sys(p,pH=pH):
            """
            Encapsulated function which iteratively sets and solves different systems
            of non-linear equations.

            Args:
                p: vector of concentrations, solved iteratively

            Returns:
                sys_eq. list of floats, solution of the equations depicted in the system.
            """

            R, T, c_H, h2o = 8.314 * 0.001 * (1 / 4.18), temp, 10 ** (- pH), 1
            sys_eq = np.zeros(len(compiled_eqs_list))
            for k,eq in enumerate(compiled_eqs_list):  # sys_eq list comprehension does not work
                sys_eq[k] = eval(eq)
            return sys_eq

        conc_i = root(non_lineal_sys, init_guess, method=solver)

        if conc_i.success == True:
            activities = Activity(conc_i.x, z_ctt, I)
            if np.isnan(activities).any() or np.isinf(activities).any():
                pass
            else:
                r2_i, rmse, mae, y_ls_values = Least_Squared(non_lineal_sys, init_guess, activities)
                if rmse < ((C_X+C_M) * threshold):
                    solved_activity_val.append(activities)
                    solved_pH_val.append(pH)
                    rmse_l.append(1/rmse)

    if len(solved_activity_val) > 0:
        solved_activity_val_T = np.array(solved_activity_val).T
    else:
        solved_activity_val_T = list()

    Kf_dft = screen_log_Kf_BiMetal(solved_activity_val_T, solved_pH_val, v_ctt, ref_idxs=ref_idx)

    return Kf_dft

def Speciation_from_Formation_singlemetal(lgkf,C, pH_grid, labels, ref_stoich, solver='hybr',acc_thr=5, temp=298.15):
    '''Computes the speciation diagram for a given speciation model from its corresponding
    formation constants (lgkf)
    Args:
        lgkf: list of floats, logarithmic formation constants log10(Kf) for all species.
        C: float, total metal concentration.
        pH_grid: array of floats, pH range to compute speciation.
        labels: list of strings, labels of the species in the diagram.
        ref_stoich: tuple of integers, stoich. coefs. for M, O and H of the reference species.
        solver: string, method to be used for nonlinear system resolution from SciPy.
        acc_thr: integer, no. of iterations used to recompute problematic pH values.
        temp: float, temperature in Kelvin
    Returns:
        solved_pH_val: list of floats, pH values.
        solved_activity_val_T: 2D NumPy array of floats, Nspecies x NpH, of concentration values.
    '''

    init_guess = [np.zeros(len(labels))]
    v_ctt = [Lab_to_stoich(lab) for lab in labels]
    M_Ratio = [v[0] for v in v_ctt]
    solved_pH_val,solved_concentration_val = list(),list()

    ind_ref = v_ctt.index(ref_stoich)

    eq1 = "p[{a}] * h2o ** 0 - (10 ** {b}) * (c_H ** {Q}) * (p[{c}] ** {P})"
    compiled_eqs_list = list()
    for ind, lg in enumerate(lgkf):
        if ind == ind_ref:
            pass
        else:
            coeff = Coefficients_Formation(v_ctt[ind_ref], v_ctt[ind])
            P, Q = coeff[0], coeff[1]
            a = ind  # Products
            c = ind_ref  # Reference
            eq = eq1.format(a=a, b=lg, c=c, Q=Q, P=P)
            compiled_eqs_list.append(compile(eq, '<string>', 'eval'))

    mass_eq = ''.join([str(M_Ratio[i]) + " * p[" + str(i) + "] + " for i in range(len(M_Ratio))]) + str(-C)
    compiled_eqs_list.append(compile(mass_eq, '<string>', 'eval'))
    for idx, pH in enumerate(pH_grid):
        def non_lineal_sys(p,pH=pH):
            """
            Encapsulated function which iteratively sets and solves different systems
            of non-lineal equations.

            Args:
                p. List of floats, concentrations of all species
            Returns:
                sys_eq. list of floats, solution of the equations depicted in the system.
            """

            R, T, c_H, h2o = 8.314 * 0.001 * (1 / 4.18), temp, 10 ** (- pH), 1
            equations_list, sys_eq = list(), list()

            for eq in compiled_eqs_list:  # sys_eq list comprehension does not work
                sys_eq.append(eval(eq))
            sys_eq = np.array(sys_eq)

            return sys_eq
        conc_i = root(non_lineal_sys, init_guess[idx], method=solver)

        acc = 0
        while conc_i.success == False:
            init2 = [uniform(0, C/10) for _ in range(len(lgkf))]  # Else, random numbers
            conc_i = root(non_lineal_sys, init2, method=solver)
            if acc > acc_thr:
                init.append(init2)
                conc_i.success = True
            elif conc_i.success and np.min(conc_i.x) < 0:
                conc_i.success = False
            acc += 1
        else:
            init_guess.append(conc_i.x)

        if conc_i.success:
            solved_concentration_val.append(conc_i.x)
            solved_pH_val.append(pH)

        if len(solved_concentration_val) > 0:
            solved_concentration_val_T = np.array(solved_concentration_val).T
        else:
            solved_concentration_val_T = list()

    return solved_pH_val, solved_concentration_val_T

def Speciation_from_Formation_bimetal(lgkf, C_X, C_M, pH_grid, labels, ref_stoich_X, ref_stoich_M, solver='hybr', acc_thr=5, temp=298.15):
    '''Computes the speciation diagram for a given speciation model from its corresponding
    formation constants (lgkf)
    Args:
        lgkf: list of floats, logarithmic formation constants log10(Kf) for all species.
        C_X, C_M: floats, total heteroatom and metal concentrations in the simulation.
        pH_grid: array of floats, pH range to compute speciation.
        labels: list of strings, labels of the species in the diagram.
        ref_stoich_X, ref_stoich_M: tuples of integers, stoich. coefs. for M, O and H of the reference species for heteroatom and metal.
        solver: string, method to be used for nonlinear system resolution from SciPy.
        acc_thr: integer, no. of iterations used to recompute problematic pH values.
        temp: float, temperature in Kelvin
    Returns:
        solved_pH_val: list of floats, pH values.
        solved_activity_val_T: 2D NumPy array of floats, Nspecies x NpH, of concentration values.
    '''
    solved_pH_val, solved_concentration_val, = list(), list()
    v_ctt = [Lab_to_stoich(lab) for lab in labels]
    X_Ratio = [v[0] for v in v_ctt]
    M_Ratio = [v[1] for v in v_ctt]
    ref = [v_ctt.index(ref_stoich_X), v_ctt.index(ref_stoich_M)]
    init_guess = [np.zeros(len(v_ctt))]

    compiled_eqs_list = list()

    eq1 = "p[{a}] * h2o ** {w} - (10 ** {b}) * (c_H ** {Q}) * (p[{c}] ** {P}) * (p[{d}] ** {R})"

    for ind,lg in enumerate(lgkf):
        if ind in ref:
            pass
        else:
            coeff = Coefficients_Formation_BiMetal(v_ctt[ref[0]], v_ctt[ref[1]], v_ctt[ind])
            P, R, Q, W = coeff
            a = ind  # Products
            c, d = ref
            eq = eq1.format(a=a, b=lg, c=c, d=d, w=W, Q=Q, P=P, R=R)
            compiled_eqs_list.append(compile(eq,'<string>', 'eval'))

    mass_eq1 = ''.join([str(X_Ratio[i]) + " * p[" + str(i) + "] + " for i in range(len(X_Ratio))]) + str(-C_X)
    mass_eq2 = ''.join([str(M_Ratio[i]) + " * p[" + str(i) + "] + " for i in range(len(M_Ratio))]) + str(-C_M)
    compiled_eqs_list.append(compile(mass_eq1, '<string>', 'eval'))
    compiled_eqs_list.append(compile(mass_eq2, '<string>', 'eval'))

    for idx, pH in enumerate(pH_grid):
        def non_lineal_sys(p,pH=pH):
            """
            Encapsulated function which iteratively sets and solves different systems
            of non-lineal equations.

            Args:
                p. List of floats, concentrations of all species
            Returns:
                sys_eq. list of floats, solution of the equations depicted in the system.
            """
            R, T, c_H, h2o = 8.314 * 0.001 * (1 / 4.18), temp,10 ** (- pH), 1
            sys_eq = list()
            for eq in compiled_eqs_list:  # sys_eq list comprehension does not work
                sys_eq.append(eval(eq))
            sys_eq = np.array(sys_eq)
            return sys_eq
        conc_i = root(non_lineal_sys, init_guess[idx], method=solver)
        acc = 0

        while conc_i.success == False:
            init2 = [uniform(0, (C_X+C_M)/10) for _ in range(len(lgkf))]  # Else, random numbers
            conc_i = root(non_lineal_sys, init2, method=solver)
            if acc >= acc_thr:
                init_guess.append(init2)
                conc_i.success = True
            elif conc_i.success and np.min(conc_i.x) < 0:
                conc_i.success = False
            acc += 1
        else:
            init_guess.append(conc_i.x)

        if conc_i.success:
            solved_concentration_val.append(conc_i.x)
            solved_pH_val.append(pH)

        if len(solved_concentration_val) > 0:
            solved_concentration_val_T = np.array(solved_concentration_val).T
        else:
            solved_concentration_val_T = list()

    return solved_pH_val, solved_concentration_val_T


def screen_log_Kf(solved_activity_val_T, solved_pH_val, stoich, ref_idx):
    """
    Generalized implementation of the private function log_Kf
    Computes formation constants for all compounds from a given reference,
    considering the stoichiometry of the process
    Args:
        solved_activity_val_T: list of floats, chemical compounds concentrations.
        solved_pH_val: list of floats, pH at which equations were solved.
        stoich: list of list of integers, chemical compounds stoichiometries.
        ref_idx: integer, index of the chemical compound used as reference.

    Returns:
        Kf_dft: list of floats, formation constants for all compounds
    """


    Kf_dft = list()

    for y_comp, sto in zip(solved_activity_val_T, stoich):
        kf_compound = list()
        for ph, hm, m in zip(solved_pH_val, y_comp, solved_activity_val_T[ref_idx]):
            kf_comp_i = log_Kf(hm, m, ph, sto, stoich[ref_idx])
            if kf_comp_i != None and not np.isnan(kf_comp_i) and not np.isinf(kf_comp_i):
                kf_compound.append(kf_comp_i)

        if len(kf_compound) > 0:
            Kf_dft.append(round(np.median(kf_compound), 4))
        else:
            Kf_dft.append(None)

    return Kf_dft

def log_Kf(HM, M, pH, prod_stoich, ref_stoich_M):
    """
    Private function which computes the formation constant of a single reaction.
    Note that the concentration of water is considered to be constant.

    For example:

    8 [MO4]2-  +   12 H+  ---->  [M8O26]4-  +    6H2O

    Kf = ([M8O26]4- * [H2O] ** 6 ) / ([MO4]2- ** 8 * [H+] ** 12)

    Args:
        HM: float, concentration of the product.
        M: float, concentration of monomer.
        pH: float, pH of the process
        prod_stoich: list of integers, (m, o, h) stoichiometric coefficients of the product
        ref_stoich_M: list of integers, (m, o, h) stoichiometric coefficients of the reference

    Returns:
        Kf: float, formation constant of the reaction
    """
    v_M, v_H3O, v_H2O = Coefficients_Formation(ref_stoich_M, prod_stoich)

    H2O = 1  # Concentration of H2O is constant
    H3O = 10 ** (- pH)

    _prod = np.log10(HM) + v_H2O * np.log10(H2O)
    _reac = v_M * np.log10(M) + v_H3O * np.log10(H3O)
    
    try:
        return _prod - _reac

    except ZeroDivisionError:
        return None

def log_Kf_BiMetal(HM, X, M, pH, prod_stoich, ref_stoich_X, ref_stoich_M):
    """
    Private function which computes the formation constant of a single reaction.
    Note that the concentration of water is considered to be constant.

    For example:

    2 [XO4]3- + 5 [MO4]2- + 10 H+ ---> [X2M5O23]6- + 5 H2O

    Kf = ([X2M5O23]6- * [H2O] ** 5 ) / ([MO4]2- ** 5 * [XO4]3- ** 2 * [H+] ** 10)

    Args:
        HM: float, concentration of the product.
        X: float, concentration of the heteroatom reference.
        M: float, concentration of metal reference.
        pH: float, pH of the process
        prod_stoich: list of integers, (x, m, o, h) stoichiometric coefficients of the product
        ref_stoich_X: list of integers, (x, m, o, h) stoichiometric coefficients of the reference heteroatom.
        ref_stoich_M: list of integers, (x, m, o, h) stoichiometric coefficients of the reference metal.

    Returns:
        Kf: float, formation constant of the reaction
    """
    v_X, v_M, v_H3O, v_H2O = Coefficients_Formation_BiMetal(ref_stoich_X, ref_stoich_M, prod_stoich)

    H2O = 1  # Concentration of H2O is constant
    H3O = 10 ** - pH

    _prod = np.log10(HM) + v_H2O * np.log10(H2O)
    _reac = v_X * np.log10(X) + v_M * np.log10(M) + v_H3O * np.log10(H3O)

    try:
        return _prod - _reac
    except ZeroDivisionError:
        return None

def screen_log_Kf_BiMetal(solved_activity_val_T, solved_pH_val, stoich, ref_idxs):
    """
    Generalized implementation of the private function log_Kf
    Computes formation constants for all compounds from a given reference,
    considering the stoichiometry of the process
    Args:
        solved_activity_val_T: list of floats, chemical compounds concentrations.
        solved_pH_val: list of floats, pH at which equations were solved.
        stoich: list of list of integers, chemical compounds stoichiometries.
        ref_idxs: integer list, index of the chemical compounds used as reference for heteroatom and metal.

    Returns:
        Kf_dft: list of floats, formation constants for all compounds
    """
    Kf_dft = list()

    for y_comp, sto in zip(solved_activity_val_T, stoich):
        kf_compound = list()
        for ph, hm, m1, m2 in zip(solved_pH_val, y_comp, solved_activity_val_T[ref_idxs[0]], solved_activity_val_T[ref_idxs[1]]):
            kf_comp_i = log_Kf_BiMetal(hm, m1, m2, ph, sto, stoich[ref_idxs[0]], stoich[ref_idxs[1]])
            if kf_comp_i != None and not np.isnan(kf_comp_i) and not np.isinf(kf_comp_i):
                kf_compound.append(kf_comp_i)

        if len(kf_compound) > 0:
            Kf_dft.append(round(np.median(kf_compound), 4))
        else:
            Kf_dft.append(None)

    return Kf_dft


def Least_Squared(func, init, y_values):
    """
    Evaluates the solution a system of non-linear equations on the basis
    of the Least Squares.

    Args:
        func: function to minimize.
        init: list of floats, initial solution for the system of NLEs.
        y_values: list of floats, activities of all compounds.

    Returns:
        r_value ** 2: float, determinant coefficient.
        rmse: float, Root mean squared error.
        mae: float, Mean average error.
        y_ls_values: list of floats, activities.

    """

    check = leastsq(func, init)
    y_ls_values = np.array(check[0]).tolist()

    slope, intercept, r_value, p_value, std_err = linregress(y_values, check[0])

    if np.isnan(check[0]).any() or np.isinf(check[0]).any():
        return np.inf, np.inf, np.inf, np.inf
    else:
        rmse = sqrt(mean_squared_error(y_values, check[0]))  #
        mae = mean_squared_error(y_values, check[0])
        return r_value ** 2, rmse, mae, y_ls_values
                                     

### Other calculation functions
def Ionic_Strenght(concentrations, charges):
    """
    Returns the ionic strength of a mixture. In the present function, 0.6M of NaClO4 was
    added according to the results reported by Cruywagen[1].

    Formula: I = 0.5 * sum(concentration_i * charge_i ** i)

    [1] Cruywagen, J. J. Protonation, Oligomerization, and Condensation Reactions
    of Vanadate(V), Molybdate(vi), and Tungstate(Vi). In Advances in Inorganic Chemistry;
    1999; Vol. 49, pp 127–182.

    Args:
        concentrations: list of floats, molar concentrations of all compounds.
        charges: list of integers, charges of all compounds.

    Returns:
        I: float, ionic strength of a mixture where 0.6M of NaClO4 was included.
    """

    S_factor = list()
    for c, z in zip(concentrations, charges):
        S_factor.append(c * ((z) ** 2))
    I = 0.5 * sum(S_factor)
    return I

def _Davies_Equation(I, charges):
    """
    Returns the correction factor to solutions with high concentrations at 25 °C. It is an empirical
    extension of Debye–Hückel theory which was first proposed in 1962.[1]

    [1] Davies, C. W. Ion Association; London: Butterworths, 1962

    Args:
        I: float,ionic strength of the simulation.
        charges: list of integers, charges of all compounds.

    Returns:
        F_factor. float, correction factor so as to transform concentrations to activities.
    """

    logF_factor = list()
    for z in charges:
        try:
            logF_factor.append(-0.5 * z ** 2 * (sqrt(I) / (1 + sqrt(I)) - 0.3 * I))
        except ValueError:
            logF_factor.append(0)  # WARNING!
    F_factor = [10 ** f for f in logF_factor]

    return F_factor

def log10_gamma(I, z):
    """
    Computes the decimal logarithm of the gamma factor as implemented in GEOCHEM:
    solution of the Davies equation.
    
    GEOMCHEM is also a modification of REDEQLE,
    made at the Department of Soil and Environmental
    Sciences of the University of California at Riverside. [1]


    Reference:

    [1] G. Sposito and S. V. Mattigod, "GEOCHEM: A
    Computer Program for the Calculation of
    Chemical Equilibria i n Soil Solutions and
    Other Natural Water Systems," Department of
    Soil and Environmental Sciences, University
    of California, Riverside, California (1979).

    Args:
        I: float, ionic strength of the simulation.
        z: list of integers, charges of all compounds.
    Returns
        float, computed log10 of gamma as implemented in the reference
    """

    A = 0.5116
    B = 0.32292 * pow(10, 8)

    if z == 0:  # Neutral Species
        B = 0.1
        return - B * I

    elif I <= 0.5:  # Davies Equation
        a = 1 / B
        B0 = 0.2 * A * (z) ** 2
        return - ((A * (z) ** 2 * np.sqrt(I)) / (1 + a * B * np.sqrt(I))) + B0 * I

    elif I > 0.5:  # Extension Equation
        a = (3 + abs(z)) * pow(10, -8)
        B0 = 0.041
        return - ((A * (z) ** 2 * np.sqrt(I)) / (1 + a * B * np.sqrt(I))) + B0 * I

def Activity(concentrations, charges, ionic_strength):
    """
    Returns the activity of chemical species.

    Args:
        concentrations: list of floats, molar concentrations of all compounds.
        charges: list of integers, charges of all compounds.
        ionic_strength: float, ionic strength of the simulation.
    Returns:
        activity: list of floats, activities of the chemical species.
    """

    I = ionic_strength
    logF_factor = list()

    for z in charges:
        try:
            logF_factor.append(log10_gamma(I, z))
        except ValueError:
            print("WARNING")
            logF_factor.append(0)  # WARNING!
    F_factor = [pow(10, f) for f in logF_factor]

    activity = np.array([conc * fz for conc, fz in zip(concentrations, F_factor)])

    return activity

def Coefficients_Formation(ref_stoich_M, prod_stoich):
    """It finds the stoichiometric coefficents for a general formation reactions
     with the following formula:

    [HxMyOz]ref + H+  --> [HxMyOz]prod + H2O

    Args:
        ref_stoich_M: (m, o, h) stoichiometric coefficients of the reference
        prod_stoich: (m, o, h) stoichiometric coefficients of the product
    Returns:
        sto_coeff: list of floats, computed coefficients for the formation reaction
    """

    H, W = [0 for _ in range(len(prod_stoich)-2)] + [0, 1], [0 for _ in range(len(prod_stoich)-2)] + [1, 2]


    A = np.array([ref_stoich_M, H, W])
    B = np.array(prod_stoich)
    sto_coeff = np.linalg.solve(A.T,B)

    return sto_coeff

def Coefficients_Formation_BiMetal(ref_stoich_X, ref_stoich_M, prod_stoich):
    """It finds the stoichiometric coefficents for a general formation reactions
     with the following formula:

    [HxMyOz]ref + H+  --> [HxMyOz]prod + H2O

    Args:
        ref_stoich_X: (x, m, o, h) stoichiometric coefficients of the heteroatom reference
        ref_stoich_M: (x, m, o, h) stoichiometric coefficients of the metal reference
        prod_stoich: (x, m, o, h) stoichiometric coefficients of the product
    Returns:
        sto_coeff: list of floats, computed coefficients for the formation reaction

    """
    H, W = [0 for _ in range(len(prod_stoich)-2)] + [0, 1], [0 for _ in range(len(prod_stoich)-2)] + [1, 2]

    A = np.array([ref_stoich_X, ref_stoich_M, H, W])
    B = np.array(prod_stoich)

    sto_coeff = np.linalg.solve(A.T,B)

    sto_coeff_abs = [abs(o) for o in sto_coeff]  # The sign means prod or reac but not useful in Kf
    return sto_coeff_abs

