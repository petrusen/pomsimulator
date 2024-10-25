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

#Simulation functions

def generate_graphs(adf_files, ref_compound, POM_system):
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
    if C_ref.shape:
        C0 = C_ref[m_idx]
    else:
        C0 = float(C_ref)
    return C0

def compute_lgkf_loop(R_idx, R_ene, R_type, mod_idx_vals, number_models,kwargs,
                      batch_size=1, cores=1):
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
    SuperArr = np.load(path_npz)["SupArray"]
    IndexArr = np.load(path_npz)["IndexArray"]
    C_ref = np.load(path_npz)["C_ref"]
    pHArr = np.load(path_npz)["pH"]
    labels = np.load(path_npz)["labels"]
    return SuperArr,IndexArr,C_ref,pHArr,labels

def models_sampling(sampling_type,number_models,sample_perc=10):
    if sampling_type == "random":

        mod_to_calc = int(number_models*sample_perc/100)
        print("Calculated speciation models number %d"%mod_to_calc)
        mod_idx_vals = random.sample(range(0, int(number_models)), mod_to_calc)  # Apply randomizer
        mod_idx_vals.sort()
    else:
        mod_idx_vals = list(range(0, number_models))
    return mod_idx_vals
#Speciation function

def speciation_diagram(idxs,lgkf_df,speciation_labels,pH,C_ref,ref_stoich):
    '''Wrapper function employed to solve all speciation models across a DataFrame.
    Args:
        idxs. Array of integers, indices of the models to be solved.
        lgkf_df. DataFrame containing log10(Kf) formation constants for all species.
        C. Float, maximum concentration value.
        speciation_labels. List of strings, labels of the species to solve speciation for.
        pH. List of floats, pH values.
        ref_stoich. Tuple of integers, stoich. coefs. for M, O and H of the reference species.
    Returns:
        concentrations_tuple. Tuple, containing list of indices (for bookkeeping) and a dictionary
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
        x_val, y_val_T = speciation_func(lgkf=lgkf,pH_grid=pH,labels=speciation_labels,**kwargs)
        concentrations_dict = dict()
        for lab, c in zip(speciation_labels, y_val_T):
            concentrations_dict[lab] = c
    concentrations_tuple = (idxs,concentrations_dict)
    return concentrations_tuple

def read_scaling_params(scaling_path):
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

def lgkf_filtering(lgkf_df, all_idxs, scaling_params, speciation_labels):
    # Select requested indices
    if all_idxs != True:
        selected_idxs = [scaling_params["best_model"]]
        batch_size = 1
        filtered_selected_idxs = [idx for idx in selected_idxs if idx in lgkf_df.index]
        lgkf_df = lgkf_df.loc[filtered_selected_idxs, :]

    # DataFrame filtering
    lgkf_df = lgkf_df[~lgkf_df.loc[:, speciation_labels].isna().any(axis=1)]
    lgkf_df = lgkf_df * scaling_params["m"] + scaling_params["b"]
    return lgkf_df

def compute_speciation_loop(lgkf_df,speciation_labels,pH,C_ref,ref_stoich,path_to_output=None,batch_size=1,cores=1,
                            show_progress=True):
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
    Q3 = lgkf_df.quantile(0.75).median()
    minmax = (lgkf_df.max() - lgkf_df.min()).median()
    maximum = lgkf_df.max().median()

    intercept = (Q3 * mlr_coefficients[0] + minmax * mlr_coefficients[1] + maximum * mlr_coefficients[2]
                 + mlr_coefficients[3])
    return intercept


def LinearScaling(path, Labels, expKf_dict, scaling_mode="best_rmse", output_scaling="regression_output.csv",
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