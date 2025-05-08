#Standard library imports
from configparser import ConfigParser
import pkg_resources as pkgr
#Local imports
from pomsimulator.modules.helper_module import load_array, get_C0
from pomsimulator.modules.stats_module import *
from pomsimulator.modules.text_module import Print_logo


def main():
    Print_logo()
    config_file = pkgr.resource_filename(__name__, "../inputs/config_PMo.pomsim")
    config = ConfigParser()
    config.read(config_file)

    output_path = pkgr.resource_filename(__name__, config["Preparation"]["output_path"])
    system = config["Preparation"]["POM_system"]

    clust_dir = output_path + "/" + config["Clustering"]["cluster_dir"]
    path_to_npz = output_path + "/" + config["Clustering"]["npz_cluster_file"]

    groups_file_path = clust_dir + "/groups.csv"
    sel_groups_idx = [int(item) for item in config["Clustering"]["sel_groups"].split(",")]
    features_file_path = output_path + "/" + config["Clustering"]["features_file"]

    SuperArr,IndexArr,C_ref,pH,labels = load_array(path_to_npz)
    ### Apply masks
    mask = np.full(SuperArr.shape[2],True)
    for ii in range(len(system.split("_"))):
        C0 = get_C0(C_ref,ii)
        mask *= mask_models(SuperArr,labels,threshold=1.1,C=C0,m_idx=ii)
    SuperArr = SuperArr[:,:,mask]
    IndexArr = IndexArr[IndexArr != -1][mask]


    Groups_dict = dict()
    with open(groups_file_path,'r') as infile:
        lines = [line.strip().split(',') for line in infile.readlines()]
        for line in lines:
            Groups_dict[int(line[0])]= [int(element) for element in line[1:]]


    group_indices = sum([Groups_dict[jj] for jj in sel_groups_idx],[])
    NewArr = SuperArr[:,:,np.array(group_indices)]
    NewIndArr = IndexArr[np.array(group_indices)]

    np.savez_compressed(clust_dir + '/NewArr.npz', SupArray=NewArr, IndexArray=NewIndArr,
                        pH=pH, labels=labels, selected_groups=sel_groups_idx,
                        original_npz=path_to_npz,C_ref=C_ref)

    if features_file_path:
        df_full = pd.read_csv(features_file_path, index_col=0)
        df_new = df_full.iloc[group_indices,:].copy()
        df_new.to_csv(clust_dir + "/new_features.csv")

    model_idxs = [IndexArr[k] for k in sorted(group_indices)]
    with open(clust_dir + "/sel_model_indices.pomsim","w") as fmod:
        fmod.write("\n".join([str(idx) for idx in model_idxs]))

    print("Selected model indices saved to ", clust_dir)
    print("Normal termination")
if __name__ == "__main__":
    main()
