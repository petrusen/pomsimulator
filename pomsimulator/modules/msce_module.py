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
def Speciation_from_Equilibrium(idx_var, e_var, type_var,mod_idx_val, idx_ctt=None,e_ctt=None,type_ctt=None,z_ctt=None,v_ctt=None,
                                ref=None, pH_grid=None, init_guess=None, I=None, C=None, threshold=None,
                                return_rmse=True, verbose=False):
    """
    Sets multi-species chemical equilibrium provided that the reactions are
    acid-base, condensations or addition reactions. The system of non-lineal
    equations are minimized using Powell's algorithm[1].

    The parameters have been fixed to the following default values:

        ·precision=5; the number of pH values at which MSCE are solved.
        ·R=0.00198; atmospheric pressure.
        ·T=298K; temperature.
        ·[H2O]=1M; concentration of water.
        ·[Mo]total=0.005; total concentration of molybdenum.
        ·rmse <= 0.005/10; convergence criterium (10% error).

    [1] Powell, M. J. D. An Efficient Method for Finding the Minimum of a Function of Several
    Variables without Calculating Derivatives. The Computer Journal. 1964, p 155.

    Args:
        idx_var: list of integers, combination of chemical reaction indexes.
        e_var: list of floats, combination of chemical reaction energies.
        type_var: list of strings, combination of chemical reaction types.
        mod_idx_val: integer, index of the model being solved
        idx_ctt, e_ctt, type_ctt: same as three first args, but for the constant reactions (acid_base)
        z_ctt: list of integers, charges of all compounds in the dataset
        v_ctt: list of lists of integers, stoichiometric indices for chemical compounds in the dataset (metal, oxygen, hydrogen)
        ref: integer, index of the reference species in compound-related lists
        pH_grid: array of floats, all pH values to solve equations for
        init_guess: array of floats, initial guess for the equation solver
        I: float, ionic strength for the simulation
        C: float, total metal concentration in the simulation
        threshold: float, controls maximum accepted rmse for a solution
        verbose: boolean controlling the verbosity of the output

    Returns:
        Kf_dft: list of floats, computed formation constants for all compounds in the reaction network

    """
    #print("Still calculating %d"%(mod_idx_val),flush=True)
    x_val, y_val, = list(), list()
    y_ls_values = list()
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
            R, T, c_H, h2o = 8.314 * 0.001 * (1 / 4.18), 298.15, 10 ** (- pH), 1
            sys_eq =  list()

            for k,eq in enumerate(compiled_eqs_list):  # sys_eq list comprehension does not work
                sys_eq.append(eval(eq))
            sys_eq = np.array(sys_eq)
            return sys_eq

        conc_i = root(non_lineal_sys, init_guess, method='hybr')

        if conc_i.success == 1:
            activities = Activity(conc_i.x, z_ctt, I)
            if np.isnan(activities).any() or np.isinf(activities).any():
                pass
            else:
                r2_i, rmse, mae, y_ls_values = Least_Squared(non_lineal_sys, init_guess, activities)
                rmse_acc = 0
                if rmse < (C * threshold):
                    y_val.append(activities)
                    x_val.append(pH)
                    rmse_acc = rmse_acc + rmse
                rmse_l.append(rmse_acc)

        if verbose == True:
            print("pH: " + str(pH) + " === Converged: " + str(conc_i.success))

    if len(y_val) > 0:
        y_val_T = [[y_val[j][i] for j in range(len(y_val))] for i in range(len(y_val[0]))]
    else:
        y_val_T = list()

    Kf_dft = screen_log_Kf(y_val_T, x_val, v_ctt, ref_idx=ref)

    return Kf_dft #, rmse_l

def screen_log_Kf(y_val_T, x_val, stoich, ref_idx):
    """
    Generalized implementation of the private function log_Kf
    Computes formation constants for all compounds from a given reference,
    considering the stoichiometry of the process
    Args:
        y_val_T: list of floats, chemical compounds concentrations.
        x_val: list of floats, pH at which equations were solved.
        stoich: list of list of integers, chemical compounds stoichiometries.
        ref_idx: integer, index of the chemical compound used as reference.

    Returns:
        Kf_dft: list of floats, formation constants for all compounds
    """


    Kf_dft = list()
    deleteme, pph = list(), list()

    for y_comp, sto in zip(y_val_T, stoich):
        kf_compound = list()
        for x, hm, m in zip(x_val, y_comp, y_val_T[ref_idx]):
            kf_comp_i = log_Kf(hm, m, x, sto, stoich[ref_idx])
            if kf_comp_i != None and not np.isnan(kf_comp_i) and not np.isinf(kf_comp_i):
                kf_compound.append(kf_comp_i)
                pph.append(x)
        if len(kf_compound) > 0:  # and not np.isinf(kf_comp_i):
            Kf_dft.append(round(sum(kf_compound) / len(kf_compound), 4))
        else:
            Kf_dft.append(None)

    return Kf_dft

def log_Kf(HM, M, pH, stoich, ref_stoich):
    """
    Private function which computes the formation constant of a single reaction.
    Note that the concentration of water is considered to be constant.

    For example:

    8 [MoO4]2-  +   12 H+  ---->  [Mo8O26]4-  +    6H2O

    Kf = ([Mo8O26]4- * [H2O] ** 6 ) / ([MoO4]2- ** 8 * [H+] ** 12)

    pKf = -log(Kf)

    Args:
        HM: float, concentration of the product.
        M: float, concentration of monomer.
        pH: float, pH of the process
        stoich: list of integers, (m, o, h) stoichiometric coefficients of the product
        ref_stoich: list of integers, (m, o, h) stoichiometric coefficients of the reference

    Returns:
        Kf: float, formation constant of the reaction
    """
    v_M, v_H3O, v_H2O = Coefficients_Formation(ref_stoich, stoich)

    H2O = 1  # Concentration of H2O is constant
    H3O = 10 ** (- pH)

    _prod = np.log10(HM) + v_H2O * np.log10(H2O)
    _reac = v_M * np.log10(M) + v_H3O * np.log10(H3O)
    
    try:
        return _prod - _reac

    except ZeroDivisionError:
        return None

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
                                     
def Speciation_from_Formation_singlemetal(C, pH, lgkf,labels,ref_stoich,solver='hybr',acc_thr=5):
    '''Computes the speciation diagram for a given speciation model from its corresponding
    formation constants (lgkf)
    Args:
        C: float, total metal concentration.
        pH: array of floats, pH range to compute speciation.
        lgkf: list of floats, logarithmic formation constants log10(Kf) for all species.
        labels: list of strings, labels of the species in the diagram.
        ref_stoich: tuple of integers, stoich. coefs. for M, O and H of the reference species.
        solver: string, method to be used for nonlinear system resolution from SciPy.
        acc_thr: integer, no. of iterations used to recompute problematic pH values.
    Returns:
        x_val: list of floats, pH values.
        y_val_T: 2D NumPy array of floats, Nspecies x NpH, of concentration values.
    '''

    init = [np.zeros(len(labels))]
    stoich = [Lab_to_stoich(lab) for lab in labels]
    M_Ratio = [v[0] for v in stoich]
    x_val,y_val = list(),list()

    ind_ref = stoich.index(ref_stoich)

    eq1 = "p[{a}] * h2o ** 0 - (10 ** {b}) * (c_H ** {Q}) * (p[{c}] ** {P})"
    compiled_eqs_list = list()
    for ind, lg in enumerate(lgkf):
        if ind == ind_ref:
            pass
        else:
            coeff = Coefficients_Formation(stoich[ind_ref], stoich[ind])
            P, Q = coeff[0], coeff[1]
            a = ind  # Products
            c = ind_ref  # Reference
            eq = eq1.format(a=a, b=lg, c=c, Q=Q, P=P)
            compiled_eqs_list.append(compile(eq, '<string>', 'eval'))

    mass_eq = ''.join([str(M_Ratio[i]) + " * p[" + str(i) + "] + " for i in range(len(M_Ratio))]) + str(-C)
    compiled_eqs_list.append(compile(mass_eq, '<string>', 'eval'))
    for idx, pHval in enumerate(pH):
        def non_lineal_sys(p):
            """
            Encapsulated function which iteratively sets and solves different systems
            of non-lineal equations.

            Args:
                p. List of floats, concentrations of all species
            Returns:
                sys_eq. list of floats, solution of the equations depicted in the system.
            """

            R, T, c_H, h2o = 8.314 * 0.001 * (1 / 4.18), 298.15, 10 ** (- pHval), 1
            equations_list, sys_eq = list(), list()

            for eq in compiled_eqs_list:  # sys_eq list comprehension does not work
                sys_eq.append(eval(eq))
            sys_eq = np.array(sys_eq)

            return sys_eq
        conc_i = root(non_lineal_sys, init[idx], method="hybr")

        acc = 0
        while conc_i.success == False:
            # print("Not converging. Step nº: " + str(acc))
            #init = [0 for _ in range(8)]
            init2 = [uniform(0, C/10) for _ in range(len(lgkf))]  # Else, random numbers

            conc_i = root(non_lineal_sys, init2, method=solver)
            acc += 1
            if acc > acc_thr:
                init.append(init2)
                break
        else:
            init.append(conc_i.x)

        if conc_i.success:
            y_val.append(conc_i.x)
            x_val.append(pH)

        if len(y_val) > 0:
            y_val_T = np.array(y_val).T
        else:
            y_val_T = list()

    return x_val, y_val_T

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

def Coefficients_Formation(ref, prod):
    """It finds the stoichiometric coefficents for a general formation reactions
     with the following formula:

    [HxMyOz]ref + H+  --> [HxMyOz]prod + H2O

    Args:
        ref: (m, o, h) stoichiometric coefficients of the reference
        prod: (m, o, h) stoichiometric coefficients of the product
    Returns:
        sto_coeff: list of floats, computed coefficients for the formation reaction
    """

    H, W = [0 for _ in range(len(prod)-2)] + [0, 1], [0 for _ in range(len(prod)-2)] + [1, 2]


    A = np.array([ref, H, W])
    B = np.array(prod)
    sto_coeff = np.linalg.solve(A.T,B)

    return sto_coeff


