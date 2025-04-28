# Standard library imports
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

# Local imports
from pomsimulator.modules.text_module import Lab_to_stoich


def mask_models(SuperArr,speciation_labels,threshold=1.1,C=0.1,m_idx=0):
    '''Retrieves a mask to exclude problematic models, with molar fractions outside
    the expected 0 - 1 range.
    Args:
        SuperArr. 3D array of floats (Nspecies x NpH x Nmodels) with all speciation information
        in concentration units.
        speciation_labels. List of strings containing labels for all species, used to retrieve
        stoichiometries.
        threshold. Float, maximum accepted molar fraction value. By default, 1.1 to avoid discarding
        values spuriously larger than unity.
        C. Float, value of total monomer concentration in mol/L. Used to go from molar fractions to
        concentrations.
        idx. Integer, position of the element in the stoichiometry used as reference. By default is 0 to select
        the metal in the MxOy-zH clusters.
    Returns:
        mask. Boolean array which is True for all models below the threshold, which can be directly used
        to index speciation arrays and feature DataFrames.
        v_ctt. List of integers containing the stoichiometric coefficient of the target species (usually, the metal)
    '''
    v_ctt = np.array([Lab_to_stoich(lab)[m_idx] for lab in speciation_labels])
    v_threshold = threshold * C / v_ctt
    mask = (SuperArr < v_threshold.reshape(-1, 1, 1)).all(axis=1).all(axis=0)
    mask2 = (SuperArr > -0.1).all(axis=1).all(axis=0)
    mask *= mask2
    return mask

### Featurization-related functions
def get_IP(arr):
    '''Computes pKa-like descriptors for peaks in the speciation diagram,
    passed as an array, detecting the inflection points of the curve.
    Args:
        arr. 2D NumPy array corresponding to a single speciation model Nspecies x NpH
    Returns:
        output. 2D NumPy array with pKa-like values for each species (Nspecies x 2) 
    '''
    last = arr.shape[1] - 1
    dx1 = np.gradient(arr,axis=1)
    dx2 = np.gradient(dx1,axis=1)
    min_d1 = np.argmin(dx1,axis=1)
    max_d1 = np.argmax(dx1,axis=1)
    output = []
    for ii,(idx1,idx2) in enumerate(zip(min_d1,max_d1)):
        if (idx1 == last):
            idx1 -= 1
        if (idx2 == last):
            idx2 -= 1
        if (idx1 == 0):
            idx1 += 1
        if (idx2 == 0):
            idx2 += 1

        p1l = dx2[ii,idx1-1]
        p1r = dx2[ii,idx1+1]
        p2l = dx2[ii,idx2-1]
        p2r = dx2[ii,idx2+1]

        if (p1l*p1r < 0):
            ip1 = idx1
        else:
            ip1 = arr.shape[1] - 1

        if (p2l*p2r < 0):
            ip2 = idx2
        else:
            ip2 = 0

        output.extend([ip1,ip2])
    return np.array(output)

def get_features(model,pH):
    '''Computes features for a speciation diagram: peak width, peak center position,
    logarithm of the peak height and area.
    Args:
        model. 2D NumPy array corresponding to a single speciation model Nspecies x NpH
        pH. 1D NumPy array with all pH values.
        v_ctt. List of integers containing the stoichiometric coefficient of the target species (usually, the metal)
    Returns:
        out_array. 1D NumPy array (4*Nspecies) containing the four features for each of the species
        in the system.
    '''
    indices = get_IP(model)
    IP_vals = pH[indices]
    width = np.abs(IP_vals[1::2] - IP_vals[0::2])
    pos = IP_vals[0::2] + width*0.5
    peaks = np.abs(pH.reshape(1,model.shape[1]) - pos.reshape(model.shape[0],1))
    peak_idx = np.argmin(peaks,axis=1)
    height = model[range(model.shape[0]),peak_idx]
    log_height = -np.log10(height)
    log_height2 = np.nan_to_num(log_height,nan=30,neginf=30,posinf=30)
    area = np.abs(np.trapz(model, x=pH, axis=1))
    out_array = np.concatenate([width,pos,log_height2,area])
    return out_array

def get_features_array(SuperArr,speciation_labels,pH):
    '''Computes features for an array of speciation diagrams
    Args:
        SuperArr. 3D array of floats (Nspecies x NpH x Nmodels) with all speciation information
        in concentration units.
        speciation_labels. List of strings containing labels for all species, used to retrieve
        stoichiometries.
        pH. 1D NumPy array with all pH values.
        v_ctt. List of integers containing the stoichiometric coefficient of the target species (usually, the metal)
    Returns:
        df. DataFrame with computed features for all models, with shape Nmodels x 4*Nspecies
    '''
    feat_list = [get_features(SuperArr[:,:,k],pH) for k in range(SuperArr.shape[2])]
    feature_names = [name + "-" + lab for name in ["width", "pos", "height","area"] for lab in speciation_labels]
    df = pd.DataFrame(feat_list, columns=feature_names)
    return df

### Clustering-related functions
def get_clusters(data, n_clusters, normalize=True):
    '''Clusterize a DataFrame containing features for a set of speciation models,
    using K-Means with predefined n_clusters and a PCA for visualization of the
    dimensionality reduction.
    Args:
        data. DataFrame with computed features for all models, with shape Nmodels x 4*Nspecies
        n_clusters. Integer, number of clusters to use for K-Means.
        normalize. Boolean, if True normalize all information in the DF using sklearn.StandardScaler
    Returns:
        cluster_dict. Dictionary containing clustering information:
            kmeans. Fitted K-Means object from scikit-learn.
            PCA: Fitted PCA object from scikit-learn.
            x_pca. 2D NumPy array containing features transformed by PCA.
    '''
    if normalize:
        scaler = StandardScaler()
        work_data = scaler.fit_transform(data)
    else:
        work_data = data
    pc = PCA()
    xpca = pc.fit_transform(work_data)
    km = KMeans(n_clusters=n_clusters)
    km.fit(work_data)
    cluster_dict = {"kmeans":km,"PCA":pc,"x_pca":xpca}
    return cluster_dict

def get_cluster_members(cluster_info,verbose=False):
    '''Get the indices of the elements pertaining to each group generated by the
    K-Means algorithm.
    Args:
        km_obj. Fitted K-Means object from scikit-learn generated by get_clusters.
        verbose. Boolean, if True show no. of elements per group through stdout.
    Returns:
        groups. List of lists of integers, containing model indices in each K-Means-based group.

    '''
    km_obj = cluster_info["kmeans"]
    tags = sorted(set(km_obj.labels_))
    groups = []
    for tgt in tags:
        sel = [idx for idx,tag in enumerate(km_obj.labels_) if tag == tgt]
        groups.append(sel)
        if verbose:
            print("group %d, %d elements" % (tgt,len(sel)))
    return groups

def write_group_file(outfile,groups):
    '''Generate a text file containing indices of the models in each group from K-Means
    Args:
        outfile. String, name of the file to save information to.
        groups. List of lists of integers, containing model indices in each K-Means-based
        group as generated by get_cluster_members.
    Returns:
        None, generates text file.
    '''
    with open(outfile, 'w') as fout:
        for ii, grp in enumerate(groups):
            fout.write(str(ii)+','+','.join([str(gr) for gr in grp]) + "\n")
            print("Cluster %d, %d elements" % (ii, len(grp)))
    return None

### Boxplot-filtering related functions
def get_boxplot_data(arr):
    '''Computes boxplot parameters (quantiles, median, IQR, whiskers and outliers)
    for a given 1D array, which in this case will correspond to a slice of size Nmodels
    of the speciation array for a single species at a single pH.
    Args:
        arr. 1D NumPy array with numeric values

    Returns:
        bp_dict. Dictionary containing boxplot parameters: quartiles, median, IQR,
    whisker limits, no. of points in the box, no. of points in the whiskers and no. of outliers
    '''
    q1, q3 = [np.percentile(arr, perc) for perc in [25, 75]]
    q2 = np.median(arr)
    iqr = q3 - q1
    w_up_val = q3 + 1.5 * iqr
    w_down_val = q1 - 1.5 * iqr
    w_up = arr[(np.abs(arr - w_up_val)).argmin()]
    w_down = arr[(np.abs(arr - w_down_val)).argmin()]
    bp_dict = {"q1": q1, "median": q2, "q3": q3, "iqr": iqr, "wup": w_up, "wdown": w_down}

    # Count outliers & points in box
    in_box = len(arr[(arr >= q1) & (arr <= q3)])
    in_whisk = len(arr[(arr >= w_down) & (arr <= w_up)])
    outliers = len(arr[(arr > w_up) | (arr < w_down)])
    bp_dict["Nbox"] = in_box
    bp_dict["Nwhisk"] = in_whisk
    bp_dict["Noutliers"] = outliers

    return bp_dict

def get_bad_models(pH_array, super_array, target_pH, target_mol_idx):
    '''Selects models corresponding to outliers for a given molecule at a given pH
    in the 3D array containing a set of speciation diagrams, and returns the
    corresponding indices so they may be discarded.
    Args:
        pH_array. 1D NumPy array with all pH values.
        super_array. 3D array of floats (Nspecies x NpH x Nmodels) with all speciation information
        in concentration units.
        target_pH. Float, pH value to be selected.
        target_mol_idx. Integer, index of the species to be selected.
    Returns:
        models_out. List of integers containing the models to be discarded (outliers)
        for the selected species and pH.
    '''
    pH_idx = np.where(abs(pH_array - target_pH) < 0.01)[0][0]
    xdata = super_array[target_mol_idx, pH_idx, :]
    bp = get_boxplot_data(xdata)
    models_out = np.argwhere((xdata > bp["wup"]) | (xdata < bp["wdown"])).reshape(-1)
    return list(models_out)

def boxplot_filtering(ndx, pH, labels, super_array):
    '''Remove all models with outliers for a given species throughout an array of
    speciation models, returning the array that only contains the models inside the
    whiskers of a boxplot (+- 1.5 times the IQR of the sample).
    Args:
        ndx. Integer, index of the species to be selected.
        pH. 1D NumPy array with pH values to be used for filtering.
        labels. List of strings containing labels for all species.
        super_array. 3D array of floats (Nspecies x NpH x Nmodels) with all speciation information
        in concentration units.
    Returns:
        sel_array. 3D array of floats (Nspecies x NpH x N_validmodels) with all speciation information
        in concentration units, after excluding all the outliers for the target species in the selected
        range of pH
    '''
    Nmodels = super_array.shape[2]
    target_spc = labels[ndx]
    # Retrieve the set of "bad" models for every pH and join the corresponding sets
    bad_model_sets = [set(get_bad_models(pH, super_array, phx, ndx)) for phx in pH]
    all_out_models = set().union(*bad_model_sets)
    valid_models = [idx for idx in range(Nmodels) if idx not in all_out_models]
    sel_array = super_array[:, :, np.r_[valid_models]]
    return sel_array
