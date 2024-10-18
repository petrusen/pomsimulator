import numpy as np
import os
from multiprocessing import Pool, cpu_count
import time
from itertools import repeat

# Local imports
from pomsimulator.modules.text_module import Print_logo,Read_csv,Lab_to_stoich,write_speciationparameters
from pomsimulator.modules.msce_module import Speciation_from_Formation_singlemetal,Speciation_from_Formation_bimetal,starmap_with_kwargs
from pomsimulator.modules.DataBase import *

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
        np.savez_compressed(path_to_output, SupArray=FilteredSuperArray, IndexArray=IndexArray, pH=pH)
    return FilteredSuperArray,IndexArray