import numpy as np
import os
from multiprocessing import Pool, cpu_count
import time
from itertools import repeat,compress,islice,product
import random
import pandas as pd
# Local imports
from pomsimulator.modules.text_module import Print_logo,Read_csv,Lab_to_stoich,write_speciationparameters,Bader_Parser
from pomsimulator.modules.msce_module import *
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.graph_module import *
from pomsimulator.modules.stats_module import mask_models,get_boxplot_data
from pomsimulator.modules.plotting_module import Reaction_Map_2D_monometal,Reaction_Map_3D_monometal

#Simulation functions

def generate_graphs(adf_files, ref_compound, POM_system):
    '''
    Generates molecular graphs for a list of ADF2019 files containing QTAIM results for connectivity
    Args:
        adf_files. list of strings, path to ADF output files
        ref_compound: string or list of strings, labels of the reference compound(s) for either a metal or a
        heteroatom and a metal.
        POM_system: string defining the POM system: element name for IPAs and X_M string for HPAs (heteroatom-atom)
    Returns:
        G1_list: list of molecular graphs generated from QTAIM connectivity
        G1_labels: list of strings, labels of the molecular set, of the form MbbOcc-dH or XaaMbbOcc-dH
        graphs_info: dictionary containing information parsed from the molecular set
            - z_ctt: list of integers, charges for all molecules in the set
            - v_ctt: list of tuples of integers, stoichiometries (bb,cc,d) or (aa,bb,cc,d) for all molecules in the set
            - water: dictionary mapping water-based molecule names to their Gibbs free energies
            - ref_idx: integer (IPA) or list of integers (HPA), containing indices of the reference species for metal or heteroatom and metal
            - num_molec: integer, number of molecules in the set
            - compounds_set: list of lists, number of molecules for each nuclearity
            - unique_labels: list of strings, names of the ADF files for each nuclearity

    '''
    G1_list, water, G1_labels = list(), dict(), list()
    for idx, f in enumerate(adf_files):
        adf_dict = Bader_Parser(f)
        label = adf_dict['label']
        if label in ['H3O', 'H2O', 'H5O2', 'H4O2']:
            water[label] = adf_dict['Gibbs']
        else:
            Gi = Molecule_to_Graph(idx, **adf_dict)
            G1_list.append(Gi)
            G1_labels.append(label)
    if isinstance(ref_compound,list):
        ref_idx = [G1_labels.index(ref) for ref in ref_compound]
    else:
        ref_idx = G1_labels.index(ref_compound)
    print("Length:", len(G1_labels), ref_idx, G1_labels)
    elements = POM_system.split("_")
    Z = [Z_dict[elem] for elem in elements]
    valence = [valence_dict[elem] for elem in elements]
    charges = Molecule_Charge(G1_list, Z, valence)
    stoich = Molecule_Stoichiometry(G1_list, Z)
    compounds_set, unique_labels = Create_Stoich(G1_labels)
    num_molec = len(G1_list)
    graphs_info = {"z_ctt": charges, "v_ctt": stoich, "water": water, "ref_idx": ref_idx, "num_molec": num_molec,
                   "compounds_set": compounds_set, "unique_labels": unique_labels}
    return G1_list, G1_labels, graphs_info

def get_C0(C_ref,m_idx):
    '''Convenience function to retrieve concentration values automatically for either IPA or HPA systems
    Args:
        C_ref: numpy array, either 0D (IPA) or 1D (HPA), containing the initial concentration(s) of metal or metal and heteroatom
        m_idx: integer. Marks the element whose concentration is selected: in XxMmOn HPA species, 0 -> X, 1 -> M.
    Returns:
        C0: float, concentration of the selected reference species
    '''
    if C_ref.shape:
        C0 = C_ref[m_idx]
    else:
        C0 = float(C_ref)
    return C0

def compute_lgkf_loop(R_idx, R_ene, R_type, mod_idx_vals, number_models, kwargs,
                      batch_size=1, cores=1):
    '''Wrapper function for the calculation of formation constants for a set of speciation models.
    Args:
        R_idx: list of lists of integers, chemical reaction indexes organized by nuclearity
        R_ene: list of lists of floats, chemical reaction energies organized by nuclearity
        R_type: list of lists of strings, chemical reaction types organized by nuclearity
        mod_idx_vals: list of integers, speciation model numbers to be solved
        kwargs: dictionary containing parameters required for speciation: for IPA Speciation_from_Equilibrium
                for HPA Speciation_from_Equilibrium_bimetal
        batch_size: integer, size of the batch of models to be solved at a time
        cores: integer, number of cores used for parallel calculation

    Returns:
        data: list of lists of logKf values for every solved speciation model, sorted according to labels' order

    '''
    if isinstance(kwargs["ref_idx"],list):
        speciation_func = Speciation_from_Equilibrium_bimetal
    else:
        speciation_func = Speciation_from_Equilibrium
    data = list()
    # Set up speciation models
    models_to_explore = set(mod_idx_vals)
    n_batches = int(len(mod_idx_vals) / batch_size)
    print("Number of batches = %d" % n_batches)
    _idx_var, _e_var, _type_var = product(*R_idx), product(*R_ene), product(*R_type)
    bool_sample = (idx in models_to_explore for idx in range(number_models))
    models_to_calculate = compress(zip(_idx_var, _e_var, _type_var), bool_sample)
    var = list()
    acc = 0
    for obj in models_to_calculate:
        var_args = (obj[0], obj[1], obj[2])
        var.append(var_args)
        acc += 1

    kwargs_iter = repeat(kwargs)

    print("Enter formation constant calculation")

    for idx in range(n_batches + 1):
        t0 = time.time()

        low_lim = batch_size * idx
        up_lim = batch_size * (idx + 1)
        if idx == n_batches:
            current_var = var[low_lim:]
        else:
            current_var = var[low_lim:up_lim]
        args_iter = current_var

        with Pool(cores) as ThreadPool:  # HPC
            data = data + starmap_with_kwargs(ThreadPool, speciation_func,
                                              args_iter, kwargs_iter)
        t1 = time.time()
        progress = idx * 100 / n_batches
        named_tuple = time.localtime()  # get struct_time
        time_string = time.strftime("%m/%d/%Y, %H:%M:%S", named_tuple)
        print(time_string,
              "[" + "".join(
                  ['#' if i < progress else " " for i in range(0, 100, 2)]) + "]" + " progress=%6.3f" % progress,
              "time of batch = %.2f s" % (t1 - t0))
    return data

def load_array(path_npz):
    '''Reads a NPZ-formatted array containing speciation information, as produced by
    Speciation_from_Formation_singlemetal or Speciation_from_Formation_bimetal
    Args:
        path_npz: string, path to the NPZ file
    Returns:
        SuperArr: 3D NumPy array of size Nspc x NpH x Nmodels containing concentration values.
        IndexArr: 1D NumPy array of size Nmodels containing model indices to map original indices to positions in SuperArr
        C_ref: 0D or 1D NumPy array containing initial concentration(s) of the reference species
        pHArr: 1D NumPy array of size NpH containing pH values used for the speciation calculation
        labels: 1D NumPy array of size Nspc containing species' labels

    '''
    SuperArr = np.load(path_npz)["SupArray"]
    IndexArr = np.load(path_npz)["IndexArray"]
    C_ref = np.load(path_npz)["C_ref"]
    pHArr = np.load(path_npz)["pH"]
    labels = np.load(path_npz)["labels"]
    return SuperArr,IndexArr,C_ref,pHArr,labels

def models_sampling(sampling_type,number_models,sample_perc=10):
    '''Samples speciation models from a total population
    Args:
        sampling_type: string, type of sampling to be used, either "random" or "all" to get all models
        number_models: integer, total number of models in the population
        sample_perc: percentage of models to be included in the random sample
    Returns:
        mod_idx_vals: list of integers, speciation model numbers to be solved
    '''
    if sampling_type == "random":
        mod_to_calc = int(number_models*sample_perc/100)
        print("Calculated speciation models number %d"%mod_to_calc)
        mod_idx_vals = random.sample(range(0, int(number_models)), mod_to_calc)  # Apply randomizer
        mod_idx_vals.sort()
    else:
        mod_idx_vals = list(range(0, number_models))
    return mod_idx_vals

def speciation_diagram(idxs,lgkf_df,speciation_labels,pH,C_ref,ref_stoich):
    '''Wrapper function employed to solve all speciation models across a DataFrame.
    Args:
        idxs: array of integers, indices of the models to be solved.
        lgkf_df: DataFrame containing log10(Kf) formation constants for all species.
        C_ref: float or list of floats, initial concentration(s) of the reference species.
        speciation_labels: list of strings, labels of the species to solve speciation for.
        pH: list of floats, pH values.
        ref_stoich: tuple of integers, stoich. coefs. for M, O and H of the reference species.
    Returns:
        concentrations_tuple: tuple, containing list of indices (for bookkeeping) and a dictionary
        of concentration arrays for each species in the system
    '''
    if isinstance(C_ref,list):
        speciation_func = Speciation_from_Formation_bimetal
        C_X,C_M = C_ref
        ref_stoich_X,ref_stoich_M = ref_stoich
        kwargs = {"C_X":C_X,"C_M":C_M,"ref_stoich_X":ref_stoich_X,"ref_stoich_M":ref_stoich_M}
    else:
        speciation_func = Speciation_from_Formation_singlemetal
        kwargs = {"C": C_ref, "ref_stoich": ref_stoich}

    max_lgkf_value = 300.0

    lgkf = lgkf_df.loc[idxs,speciation_labels].to_numpy()
    if max(lgkf) > max_lgkf_value:
        concentrations_dict = []
    else:
        solved_pH_val, solved_conc_val_T = speciation_func(lgkf=lgkf,pH_grid=pH,labels=speciation_labels,**kwargs)
        concentrations_dict = dict()
        for lab, c in zip(speciation_labels, solved_conc_val_T):
            concentrations_dict[lab] = c
    concentrations_tuple = (idxs,concentrations_dict)
    return concentrations_tuple

def read_scaling_params(scaling_path):
    '''Convenience function to read the scaling parameter files generated by scale_constants.py
    Args:
        scaling_path: string, path to the file containing scale parameters
    Returns
        dictionary, containing:
            - m: float, slope of the regression
            - b: float, intercept of the regression
            - best_model: integer, index of the best model (if not applicable, -1)
            - mode: string, scaling mode used in scale_constants.py -> best_rmse, average, medians or universal
    '''
    with open(scaling_path, "r") as fpars:
        data_raw = [line.strip().split() for line in fpars.readlines()]
        data = data_raw[1]
        mode = data[0]
        slope = float(data[1])
        intercept = float(data[2])
        best_model = int(data[3])
        if best_model < 0:
            best_model = None
    return {"m":slope,"b":intercept,"best_model":best_model,"mode":mode}

def apply_lgkf_scaling(lgkf_df, scaling_params, speciation_labels):
    '''Filters and scales a DataFrame containing formation constants, applying the scaling parameters contained in a
       dictionary, selecting a best model if applicable and selecting requested labels
    Args:
        lgkf_df: DataFrame containing log(Kf) values produced by a simulation
        scaling_params: dictionary of scaling parameters read by read_scaling_params()
        speciation_labels: list of strings, labels used for the speciation

    Returns:
        lgkf_df: DataFrame containing scaled, selected log(Kf) values
    '''
    # Select requested indices
    if scaling_params["mode"] == "best_rmse":
        selected_idxs = [scaling_params["best_model"]]
        filtered_selected_idxs = [idx for idx in selected_idxs if idx in lgkf_df.index]
        lgkf_df = lgkf_df.loc[filtered_selected_idxs, :]
    # DataFrame filtering
    lgkf_df = lgkf_df[~lgkf_df.loc[:, speciation_labels].isna().any(axis=1)]
    lgkf_df = lgkf_df * scaling_params["m"] + scaling_params["b"]
    return lgkf_df

def compute_speciation_loop(lgkf_df,speciation_labels,pH,C_ref,ref_stoich,path_to_output=None,
                            batch_size=1,cores=1,show_progress=True):
    '''Wrapper function for the calculation of speciation diagrams for a set of speciation models.
    Args:
        lgkf_df: DataFrame containing scaled, selected log(Kf) values
        speciation_labels: list of strings, labels used for the speciation
        pH: list of floats, pH values.
        C_ref: float or list of floats, initial concentration(s) of the reference species.
        ref_stoich: list of integers (IPA) or list of lists of integers (HPA), stoich. coefs. for (X), M, O and H of the reference species.
        path_to_output: string, name of the NPZ file to store the speciation. If None, no file is created.
        batch_size: integer, size of the batch of models to be solved at a time
        cores: integer, number of cores used for parallel calculation
        show_progress: boolean, if True, print a progress bar throughout the calculation
    Returns:
        FilteredSuperArray: 3D NumPy array of size Nspc x NpH x Nmodels containing concentration values.
        FilteredIndexArray: 1D NumPy array of size Nmodels containing model indices to map original indices to positions in FilteredSuperArray

    '''
    if isinstance(C_ref,list):
        C_X,C_M = C_ref
        ref_stoich_X,ref_stoich_M = ref_stoich
        add_kwargs = {"C_ref":[C_X,C_M],"ref_stoich":[ref_stoich_X,ref_stoich_M]}
    else:
        add_kwargs = {"C_ref": C_ref, "ref_stoich": ref_stoich}

    filt_sel_idxs = list(lgkf_df.index)
    Nidx = len(filt_sel_idxs)
    ### Array initialization
    SuperArray = np.zeros((len(speciation_labels), len(pH), Nidx))
    IndexArray = np.full((len(filt_sel_idxs)), -1, dtype=int)
    k = 0

    # Batch setup and argument preparation
    n_batches = int(Nidx / batch_size)
    kwargs = dict(lgkf_df=lgkf_df,speciation_labels=speciation_labels, pH=pH)
    kwargs.update(add_kwargs)
    print("Number of batches = %d" % n_batches)

    # Main calculation loop
    for idx in range(n_batches + 1):
        t0 = time.time()
        list_concent = list()
        if idx == n_batches:
            batch = filt_sel_idxs[batch_size * idx:]
        else:
            batch = filt_sel_idxs[batch_size * idx:batch_size * (idx + 1)]

        args_iter = [[item] for item in batch]
        kwargs_iter = repeat(kwargs)
        with Pool(cores) as ThreadPool:  # HPC
            list_concent = list_concent + starmap_with_kwargs(ThreadPool, speciation_diagram, args_iter,
                                                              kwargs_iter)

        for idxs, concent in list_concent:
            if len(concent) > 0:
                if len(concent[speciation_labels[0]]) != len(pH):
                    print("Shape concent is different from shape pH for idx:", idxs)
                else:
                    for j, (molec, conc) in enumerate(concent.items()):
                        SuperArray[j, :, k] = conc
                    IndexArray[k] = idxs
                    k += 1
        t1 = time.time()
        progress = idx * 100 / n_batches
        named_tuple = time.localtime()  # get struct_time
        time_string = time.strftime("%m/%d/%Y, %H:%M:%S", named_tuple)
        if show_progress:
            print(time_string, "[" + "".join(
                ['#' if i <= progress else " " for i in range(0, 100, 2)]) + "]" + " progress=%5.2f" % progress,
                  "time of batch = %.3f s" % (t1 - t0))

    FilteredSuperArray = SuperArray[:, :, 0:k].astype(np.float32)
    if path_to_output:
        np.savez_compressed(path_to_output, SupArray=FilteredSuperArray, IndexArray=IndexArray,C_ref=C_ref,
                            pH=pH,labels=speciation_labels)
    return FilteredSuperArray,IndexArray

def Internal_Lab_Gen(df, lgkf_dict):
    '''Determines the species having experimental data from the dictionary
    to filter out the DataFrame.
    Args:
        df. Pandas DataFrame containing lgkf values for a set of speciation models.
        lgkf_dict. Dict mapping species' labels to experimental rate constants.
    Returns:
        filt_df. Filtered DataFrame only having lgkf values for species with experimental values.
        int_labels. List of strings, labels of all species with experimental constants.
        Kexp. NumPy array of floats, exp. constant values.

    '''
    int_labels = list(lgkf_dict.keys())
    filt_df = df.loc[:, int_labels]

    filt_df = filt_df[~filt_df.isna().any(axis=1)]

    Kexp = [lgkf_dict[lab] for lab in int_labels]
    Kexp = np.array(Kexp)

    return filt_df, int_labels, Kexp

def lgKf_scaling(filt_df, Kexp):
    '''Applies linear regression to the whole set of speciation models in the DataFrame against the experimental
    constants, returning a DF of regression parameters for each model.
    Args:
        filt_df. Pandas DataFrame as generated by Internal_Lab_Gen, only having columns for species with exp. constant.
        Kexp. NumPy array of floats, exp. constant values.
    Returns:
        df2. Pandas DataFrame with slope, intercept, regression coefficient, standard error and rmse for each model.
    '''
    result = filt_df.apply(func=linregress, axis=1, y=Kexp)
    df2 = pd.DataFrame.from_dict(dict(zip(result.index, result.values)), orient="index",
                                 columns=["m", "b", "r", "p", "std_err"])
    a2 = df2.to_numpy()
    arr = filt_df.to_numpy()
    arr_sc = arr * a2[:, 0].reshape(-1, 1) + a2[:, 1].reshape(-1, 1)
    rmse = np.sqrt(np.apply_along_axis(func1d=mean_squared_error, axis=1, arr=arr_sc, y_pred=Kexp))
    df2.loc[:, 'rmse'] = rmse
    return df2

def df_2_boxplot(df):
    '''Complete workflow to generate boxplot values for computed DFT constants in a dataset
   and plot them in the positions marked by experimental constants
    Args:
        df. Pandas DataFrame containing lgkf values for a set of speciation models.
        Exp_kf_dict. Dict mapping species' labels to experimental rate constants.
        color. String, color specification for matplotlib.
        box_height. Float, height of boxplot for plotting.
        ax. Matplotlib Axis object to put figure in. If None, generate new figure.
        remove_outliers. Boolean, if True, filter out all the outliers in the boxplot.
    Returns:
        boxplot_dict. Dict indexed by labels of species with all information for each
        boxplot as produced by get_boxplot_data
        Produces figure with all boxplots.
    '''
    species = df.columns
    boxplot_dict = {}
    for ii, sp in enumerate(species):
        arr = df.loc[:, sp].to_numpy()
        boxplot = get_boxplot_data(arr)
        boxplot_dict[sp] = boxplot
    return boxplot_dict

def apply_MLR(lgkf_df):
    '''Apply MLR parameters (available in DataBase) derived from universal POM scaling to a set of log(Kf) values
    Args:
        lgkf_df: DataFrame containing log(Kf) values produced by a simulation
    Returns:
        intercept: float, predicted value of the intercept for the set of DFT constants
    '''
    Q3 = lgkf_df.quantile(0.75).median()
    minmax = (lgkf_df.max() - lgkf_df.min()).median()
    maximum = lgkf_df.max().median()

    intercept = (Q3 * mlr_coefficients[0] + minmax * mlr_coefficients[1] + maximum * mlr_coefficients[2]
                 + mlr_coefficients[3])
    return intercept

def LinearScaling(path, expKf_dict, scaling_mode="best_rmse", output_scaling="regression_output.csv",
                  Metal=None, output_path="."):
    '''Wrapper function to apply linear scaling to a set of lgkf values in a
    CSV file.
    Args:
        path. String, path to the file to read constants from.
        Labels. List of strings, complete list of labels for the system being treated.
        d_Cruy. Dict mapping species' labels to experimental rate constants.
        scaling_mode. String, can be:
            - best_rmse. Gather the model with the lowest RMSE and select its m, b parameters.
            - average. Get the average m and b parameters among all models.
        universal_scaling. Boolean, if True consider boxplots for all regressions, else
        only plot best models.
    Returns:
        None: generates image and output files.
            - regression_output.csv. CSV file with the parameters of all individual regressions.
            - scaling_params_PMo.pomsim. File containing the slope and intercept to be used in speciation.
    '''
    lgkf_df = Read_csv(path)
    if scaling_mode == "universal":
        _lgkf_df = lgkf_df
    else:
        _lgkf_df, int_labels2, Kexp = Internal_Lab_Gen(lgkf_df, expKf_dict)

    # Flag to perform all regressions
    all_regressions = scaling_mode in ["best_rmse", "average"]

    if all_regressions:
        scaling_params_df = lgKf_scaling(_lgkf_df, Kexp)
        # Sort the array and save it
        scaling_params_sorted = scaling_params_df.sort_values(by="rmse")
        scaling_params_sorted.to_csv(output_scaling)
        best_model_idx = int(scaling_params_sorted.index[0])

    if scaling_mode == "average":
        best_model_idx = -1
        rmse_average = scaling_params_df["rmse"].mean()
        slope_average = scaling_params_df["m"].mean()
        intercept_average = scaling_params_df["b"].mean()

        # Printing average of the list
        print("RMSE average =", round(rmse_average, 2), " Minimum RMSE =", round(scaling_params_df["rmse"].min(), 2))
        print("Slope average =", round(slope_average, 2))
        print("Intercept average =", round(intercept_average, 2))

        scaling_params_dict = {"m": slope_average, "b": intercept_average, "r2": None}

    elif scaling_mode == "best_rmse":
        best = scaling_params_sorted.iloc[0, :]
        scaling_params_dict = {"m": best.m, "b": best.b, "r2": best.r ** 2}
    elif scaling_mode == "medians":
        best_model_idx = -1
        boxplot = df_2_boxplot(_lgkf_df)
        Kf_dft_median = [boxplot[k]["median"] for k in int_labels2]
        Kf_exp = [expKf_dict[k] for k in int_labels2]
        reg_param = linregress(Kf_dft_median, Kf_exp)

        m = reg_param.slope
        b0 = reg_param.intercept
        r = reg_param.rvalue
        scaling_params_dict = {"m": m, "b": b0, "r2": r ** 2}

    elif scaling_mode == "universal":
        m = universal_slope
        """Apply MLR method described in xxxxxxxx"""
        b0 = apply_MLR(_lgkf_df)
        best_model_idx = -1
        scaling_params_dict = {"m": m, "b": b0, "r2": None}

    # Write to file
    if Metal == None:
        outfile = output_path + "/scaling_params.pomsim"
    else:
        outfile = output_path + "/scaling_params_%s.pomsim" % Metal
    with open(outfile, "w") as fout:
        header = "%10s %6s %6s %12s\n" % ("Mode", "m", "b", "best")
        vals = "%10s %6.4f %6.4f %12d \n" % (scaling_mode, scaling_params_dict["m"],
                                             scaling_params_dict["b"], best_model_idx)

        fout.write(header + vals)

    return scaling_params_dict

### Functions to handle nuclearities
def stoich_to_lab(at1,at2,sto):
    '''Helper function to produce strings for nuclearity stoichiometries, given a (xx,mm,oo) tuple of integers
    Args:
        at1,at2: strings, symbols of the heteroatom and atom in the target system.
        sto: tuple of integers, with stoich. coefficients of heteroatom, atom and oxygen
    Returns:
        lab: string of the form XxxMmmOoo for the target nuclearity
    '''
    lab = "%s%02d%s%02dO%02d" % (at1,sto[0],at2,sto[1],sto[2])
    return lab

def collapse_color_dict(col_dict):
    '''Simplifies a color dictionary to one without protonation states, assigning the first occurrence of that
    nuclearity in the original dict as the overall color
    Args:
        col_dict: dictionary mapping labels to color specifications.
    Returns:
        nw_col_dict: dictionary mapping nuclearity labels to color specifications.
    '''
    # remove protonation states from the color dict
    nw_col_dict = {}
    for mol,col in col_dict.items():
        nuc = re.sub("-[0-9]{1,2}H","",mol)
        if nuc not in nw_col_dict.keys():
            nw_col_dict[nuc] = col
    return nw_col_dict

def collapse_labels(labels):
    '''Simplifies a list of labels to one without protonation states, preserving the original order
    Args:
        labels: list of strings, labels of the species in the diagram.
    Returns:
        nw_labels: list of strings, labels of the nuclearities in the diagram.
    '''

    nw_labels = []
    for lab in labels:
        lab = re.sub("-[0-9]{1,2}H","",lab)
        if lab not in nw_labels:
            nw_labels.append(lab)
    return nw_labels

def nuclearity_collapser(SuperArr,Labels):
    '''Processes a Nspc x NpH x Nmodels array to join the concentrations of the protonation states of individual nuclearities, obtaining a
    Nnuc x NpH x Nmodels array
    Args:
        SuperArr: 3D NumPy array of size Nspc x NpH x Nmodels containing concentration values.
        Labels: list of strings, labels of the species in the diagram.
    Returns:
        nw_arr: array of floats, size Nnuc x NpH x Nmodels, from the addition of concentrations of the same nuclearity.
        known_nuc: list of lists of integers, containing stoichiometric indices (xx,mm,oo) for all nuclearities found in the system
    '''
    # generate stoichiometries and remove hydrogens
    stoichs = [Lab_to_stoich(lab) for lab in Labels]
    sto_noH = [tuple(sto[0:-1]) for sto in stoichs]

    # get unique values, preserving the order
    selection = defaultdict(list)
    known_nuc = []
    for ii,sto in enumerate(sto_noH):
        selection[sto].append(ii)
        if sto not in known_nuc:
            known_nuc.append(sto)
    Nnuc = len(selection.keys())

    # sum values in the array
    nw_arr = np.zeros((Nnuc,SuperArr.shape[1],SuperArr.shape[2]))
    for ii,nuc in enumerate(known_nuc):
        elems = selection[nuc]
        nw_arr[ii,:,:] = SuperArr[tuple(elems),:,:].sum(axis=0)
    return nw_arr,known_nuc

# Phase diagrams
def phase_diagram_HPA(npz_paths,v_ctt):
    """
        Gathers all concentration arrays, and selects the most predominant species at each pH value.
    Args:
        npz_paths: list of strings, paths to the NPZ files containing concentration arrays, pH, initial concentration(s) and labels
        v_ctt: List of stoichiometries for each molecule

    Returns:
        phase_diagram_X, phase_diagram_M: 2D NumPy array of shape Nrat x NpH containing the indices of the species with most % at each pair of C,pH,
                                          considering the heteroatom (X) and the metal (M)
        Ratio_list: list of floats, length Nrat, X/M ratio values for the Y-axis of the phase diagram
        pH: 1D NumPy array of shape NpH, pH values for the X-axis of the phase diagram
    """
    v_ctt_arr_X = np.array([item[0] for item in v_ctt]).reshape(-1, 1)
    v_ctt_arr_M = np.array([item[1] for item in v_ctt]).reshape(-1, 1)

    Ratio_list = list()
    for ii,path in enumerate(npz_paths):
        phase_array_dict = np.load(path)
        labels = phase_array_dict["labels"]
        C_X = phase_array_dict["C_X"]
        C_M = phase_array_dict["C_M"]
        Ratio_list.append(C_M/C_X)
        pH = phase_array_dict["pH"]
        speciation_array = phase_array_dict["SupArray"]
        if ii == 0:
            phase_diagram_X = np.zeros((len(npz_paths), len(pH)), dtype=int)
            phase_diagram_M = np.zeros((len(npz_paths), len(pH)), dtype=int)

        mask_X = mask_models(speciation_array,labels,1.1,C_X,m_idx=0)
        mask_M = mask_models(speciation_array,labels,1.1,C_M,m_idx=1)
        speciation_array_X = speciation_array[:,:,mask_X]
        speciation_array_M = speciation_array[:,:,mask_M]
        Means_X = np.mean(speciation_array_X, axis=2)*v_ctt_arr_X
        Means_M = np.mean(speciation_array_M, axis=2)*v_ctt_arr_M

        phase_diagram_X[ii, :] = np.argmax(Means_X, axis=0)
        phase_diagram_M[ii, :] = np.argmax(Means_M, axis=0)
    return phase_diagram_X,phase_diagram_M,Ratio_list,pH

def phase_diagram_IPA(npz_paths,v_ctt):
    """
        Gathers all concentration arrays, and selects the most predominant species at each pH value.
    Args:
        npz_paths: list of strings, paths to the NPZ files containing concentration arrays, pH, initial concentration(s) and labels
        v_ctt: List of stoichiometries for each molecule

    Returns:
        phase_diagram: 2D NumPy array of shape Nconc x NpH containing the indices of the species with most % at each pair of C,pH.
        C_list: list of floats, length Nconc, concentration values for the Y-axis of the phase diagram
        pH: 1D NumPy array of shape NpH, pH values for the X-axis of the phase diagram
    """
    v_ctt_arr = np.array([item[0] for item in v_ctt]).reshape(-1, 1)
    C_list = list()
    for ii,path in enumerate(npz_paths):
        phase_array_dict = np.load(path)
        C = phase_array_dict["C"]
        C_list.append(C)
        pH = phase_array_dict["pH"]
        speciation_array = phase_array_dict["SupArray"]
        if ii == 0:
            phase_diagram = np.zeros((len(npz_paths), len(pH)), dtype=int)

        mask = mask_models(speciation_array,phase_array_dict["labels"],1.1,C)
        speciation_array = speciation_array[:,:,mask]
        Means = np.mean(speciation_array, axis=2)
        Means = Means*v_ctt_arr
        phase_ratio = np.argmax(Means, axis=0)
        phase_diagram[ii, :] = phase_ratio
    return phase_diagram,C_list,pH

#CRN
def generate_CRN(G_list,stoich,reac_idx,reac_ener,reac_type,sorted_reac_idx,sorted_reac_ener,sorted_reac_type,
                 plotting_dict,plot_dict_details):
    """acid base reactions"""

    if plotting_dict["full"]:
        if plotting_dict["dimension_3d"]:
            fig,axd = Reaction_Map_3D_monometal(G_list,reac_idx,reac_ener,reac_type,stoich,
                                                All_models=True,ploting_details_dict=plot_dict_details)
            return fig,axd
        else:
            fig, ax = Reaction_Map_2D_monometal(G_list, reac_idx, reac_ener, reac_type, stoich,
                                                 All_models=True, ploting_details_dict=plot_dict_details)
            return fig,ax
    else:
        idx_ctt = list(reac_idx[0])
        e_ctt = list(reac_ener[0])
        type_ctt = list(reac_type[0])

        acc = 0
        for idx_var, e_var, type_var in zip(product(*sorted_reac_idx),product(*sorted_reac_ener),product(*sorted_reac_type)):
            if acc == plotting_dict["mod_idx"]:
                if plotting_dict["dimension_3d"]:
                    fig, axd = Reaction_Map_3D_monometal(G_list, list(idx_var) + idx_ctt, list(e_var) + e_ctt,
                                                         list(type_var)+type_ctt, stoich,
                                                         All_models=False, ploting_details_dict=plot_dict_details)
                    return fig, axd
                else:
                    fig, ax = Reaction_Map_2D_monometal(G_list, list(idx_var) + idx_ctt, list(e_var) + e_ctt,
                                                         list(type_var)+type_ctt, stoich,
                                                         All_models=False, ploting_details_dict=plot_dict_details)
                    return fig, ax
            acc += 1