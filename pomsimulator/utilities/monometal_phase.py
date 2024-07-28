# Standard library imports
import time
import numpy as np

# Local imports
from pomsimulator.modules.text_module import Read_csv,Print_logo
from pomsimulator.modules.graph_module import *
from pomsimulator.modules.msce_module import *
from pomsimulator.modules.DataBase import *

def speciationdiagram_monometal(idxs,lgkf_df,C,speciation_labels,pH,ref_stoich,stoich_lab_dict,v_ctt):
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
    max_lgkf_value = 300.0

    # print("Starting resolution of speciation model:",idxs)

    lgkf = lgkf_df.loc[idxs,speciation_labels].to_numpy()
    if max(lgkf) > max_lgkf_value:
        # print("overflow for this idx:", idxs)
        concentrations_dict = []

    else:
        x_val, y_val_T = Speciation_from_Formation_singlemetal(C, pH,lgkf,speciation_labels,ref_stoich,solver='hybr')
        concentrations_dict = dict()
        for vk, c in zip(v_ctt, y_val_T):
            concentrations_dict[stoich_lab_dict[tuple(vk)]] = c

    concentrations_tuple = (idxs,concentrations_dict)
    return concentrations_tuple



def main():
    ######################### LABELS ##############################################
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'
    Print_logo()

    # Taken from DataBase
    formation_labels = Labels_W
    speciation_labels = Labels_W_good
    ref_stoich = [1,4,0]        # WO4(2-)
    pH = np.arange(9,0,-0.01)

    C_list = np.logspace(-4,0,200)
    #C_list = np.linspace(1e-4,1,200)
    v_ctt = [Lab_to_stoich(lab) for lab in speciation_labels]
    path = "../outputs/logkf_W.txt"
    batch_size = 1
    cores = 10

    with open("../outputs/scaling_params.pomsim","r") as fpars:
        data_raw = [line.strip().split() for line in fpars.readlines()]
        data = data_raw[1]
        mode = data[0]
        slope = float(data[1])
        intercept = float(data[2])
        best_model = int(data[3])
        if best_model < 0:
            best_model = None

    selected_idxs =[]

    all_idxs = False

    ###############################################################################################
    # Read constants and scale them
    stoich_lab_dict = {tuple(Lab_to_stoich(lab)): lab for lab in formation_labels}
    _lgkfdft_df = Read_csv(path, formation_labels)

    # Subset of values only
    if all_idxs != True:
        selected_idxs = [best_model]
        batch_size = 1
        filtered_selected_lgkfdft_idxs = [idx for idx in selected_idxs if idx in _lgkfdft_df.index]
        _lgkfdft_df = _lgkfdft_df.loc[filtered_selected_lgkfdft_idxs,:]

    _lgkfdft_df = _lgkfdft_df[~_lgkfdft_df.loc[:,speciation_labels].isna().any(axis=1)]
    _lgkfdft_df = _lgkfdft_df * slope + intercept
    filtered_selected_lgkfdft_idxs = list(_lgkfdft_df.index)

    array_dict = {"conc_array":C_list,"pH":pH}

    for ii,C in enumerate(C_list):
        print("Speciation for concentration = %.6f" % C)
        SuperArray = np.zeros((len(speciation_labels), len(pH), len(filtered_selected_lgkfdft_idxs)),dtype=np.float32)
        IndexArray = np.full((len(filtered_selected_lgkfdft_idxs)),-1,dtype=int)
        t00 = time.time()
        k = 0
        n_batches = int(len(filtered_selected_lgkfdft_idxs)/batch_size)

        kwargs = dict(lgkf_df=_lgkfdft_df,C=C,speciation_labels=speciation_labels,pH=pH,ref_stoich=ref_stoich,stoich_lab_dict=stoich_lab_dict,v_ctt=v_ctt)
        print("Number of batches = %d" % n_batches)
        print("Cores = %d" % cores)
        for idx in range(n_batches + 1):
            t0 = time.time()
            list_concent = list()
            if idx == n_batches:
                batch = filtered_selected_lgkfdft_idxs[batch_size*idx:]
            else:
                batch = filtered_selected_lgkfdft_idxs[batch_size*idx:batch_size*(idx+1)]
            args_iter = [[item] for item in batch]
            kwargs_iter = repeat(kwargs)
            with Pool(cores) as ThreadPool:  # HPC
                list_concent = list_concent + starmap_with_kwargs(ThreadPool, speciationdiagram_monometal, args_iter, kwargs_iter)

            print("% solved models = ",idx*100/n_batches)

            for idxs,concent in list_concent:
                if len(concent) > 0:
                    if len(concent[formation_labels[0]]) != len(pH):
                        print("Shape concent is different from shape pH for idx:", idxs)
                    else:
                        for j, (molec, conc) in enumerate(concent.items()):
                            SuperArray[j, :, k] = conc
                        IndexArray[k] = idxs
                        k += 1

            t1 = time.time()
            print("time of batch = %.4f s" % (t1-t0))

        speciation_array = SuperArray[:,:,0:k].astype(np.float32)
        array_dict["array_%d" % ii] = speciation_array


        t11 = time.time()
        print("t = %.2f s" % (t11-t00))

    np.savez_compressed("../outputs/phase_array.npz",**array_dict)

if __name__ == '__main__':
    main()
