# Standard library imports
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns


sns.set_context('talk', font_scale=0.8)
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np

#Local imports
from pomsimulator.modules.text_module import Lab_to_stoich,Lab_to_Formula
from pomsimulator.modules.stats_module import get_boxplot_data


def Reaction_Map_3D_monometal(G1_list, all_reac_idx, all_reac_e, all_reac_type, stoich, All_models=True,ploting_details_dict=None):
    """
    Converts a graph into a reaction map. It weights the edges as a function of the
     reaction energy and uses colors to differentiate them. July 2019 Version.
     Args:

     Kwargs:

     Return:

     """

    print("Entering reaction map")
    G2_obj = nx.Graph()
    colormap_name = ploting_details_dict['colormap']
    colormap = getattr(plt.cm,colormap_name)

    trivial_yaxis, trivial_xaxis, trivial_zaxis = list(), list(), list() #These lists will contain the x,y,z positions of each graph

    """Next step will iterate over al the compounds in G1_list and according to its stoichiometry it will then give a xyz coordinates 
    (i.e. for the molybdate anion [0,1,4,0] the coordinates would be as follows: z coordinate would be the number of molybdenum atoms also called nuclearity,
    the y coordinate would be the number of protons and the x coordinate would be the ratio between oxygen atoms and molybdenum atoms plus the number of phosphorus.
    x coordinate is normalised for the trimer P00Mo03O09-0H(0,0,3))"""

    for i in range(len(G1_list)):
        list_Z = list(nx.get_node_attributes(G1_list[i], 'Z').values())
        G2_obj.add_node(i)
        if len(stoich[i]) == 3:
            trivial_zaxis.append(stoich[i][0])
            trivial_yaxis.append(stoich[i][2])
            trivial_xaxis.append(0 + (stoich[i][0]/stoich[i][1]))
        elif len(stoich[i]) == 4:
            trivial_zaxis.append(stoich[i][1])
            trivial_yaxis.append(stoich[i][3])
            trivial_xaxis.append(stoich[i][0] + (stoich[i][1]/stoich[i][2]))
    if All_models == True: #when functions is called, arguments can be all reactions or the reactions for 1 particular model.
        reac_idx = [item for sublist in all_reac_idx for item in sublist]
        reac_e = [item for sublist in all_reac_e for item in sublist]
        reac_type = [item for sublist in all_reac_type for item in sublist]
    else:
        reac_idx = all_reac_idx
        reac_e = all_reac_e
        reac_type = all_reac_type

    '''When edges are created, they are given an attribute that represents the reaction energy of that edge.
        It will alter be used to colour the reaction map depending on energies'''

    for i in range(len(reac_idx)):
        G2_obj.add_edge(reac_idx[i][0], reac_idx[i][-1], attri=reac_e[i])
        if len(reac_idx[i]) > 2:
            G2_obj.add_edge(reac_idx[i][0], reac_idx[i][-2], attri=reac_e[i])

    nodes = nx.nodes(G2_obj)
    edges = nx.edges(G2_obj)

    coloredges = list()
    for i, j in edges:
        coloredges.append(G2_obj[i][j]['attri'])

    mosaic = """AB
                AB"""
    plt.rcParams['grid.color'] = "#e5e5e5"


    fig, axdict = plt.subplot_mosaic(mosaic,figsize=(6.4*2.3,5.5*2.3),gridspec_kw={"width_ratios":[0.1,1],
                                                                                  "wspace":-0.2})
    b_spec = axdict['B'].get_subplotspec()
    axdict['B'].remove()
    axdict['B'] = fig.add_subplot(b_spec, projection='3d')

    pos_l = list()

    for idx, n in enumerate(nodes):
        pos_l.append((n, (trivial_xaxis[idx], trivial_yaxis[idx], trivial_zaxis[idx])))
    pos = dict(pos_l)
    pos_edges = np.array([(pos[u],pos[v]) for u,v in G2_obj.edges()])
    pos_nodes = np.array([item[1] for item in pos_l])

    for i in range(len(pos_nodes)):
        axdict["B"].scatter(*pos_nodes.T,color=ploting_details_dict['node_color'],s=50)
        norm = matplotlib.colors.Normalize(vmin=min(coloredges), vmax=max(coloredges))
        edgecolor = [colormap(norm(e_edge)) for e_edge in coloredges]
    for i,edge in enumerate(pos_edges):
        axdict["B"].plot(*edge.T,color=edgecolor[i])

    col_bar = axdict["B"].scatter(coloredges, coloredges, s=0, cmap=colormap, c=coloredges)
    plt.colorbar(col_bar, ax=axdict["A"],label='Reaction Energy',location='left')

    axdict["B"].set_xlim([min(trivial_xaxis),max(trivial_xaxis)])
    axdict["B"].set_ylim([min(trivial_yaxis),max(trivial_yaxis)])
    axdict["B"].set_xlabel(ploting_details_dict['x_axis_lab'], fontsize=15,labelpad=20)
    axdict["B"].set_ylabel(ploting_details_dict['y_axis_lab'], fontsize=15,labelpad=20)
    axdict["B"].set_zlabel(ploting_details_dict['z_axis_lab'], fontsize=15,labelpad=20)
    axdict["B"].set_title(ploting_details_dict['plot_title'], fontsize=16,y=0.9)

    axdict['A'].axis("off")



        # for ii in np.arange(270, 630, 1):
        #     axdict["B"].view_init(elev=0, azim=ii)
        #     plt.savefig("/home/jbuils/kimikhome/POMSimulator/KEGs_simulator/examples/figures/3d_map_full_rotate_cyberpunk2/3dmap%d.png" % ii,transparent=False)

    return fig,axdict

def Reaction_Map_2D_monometal(G1_list, all_reac_idx, all_reac_e, all_reac_type, stoich, All_models=True,ploting_details_dict=None):
    """
    Converts a graph into a reaction map. It weights the edges as a function of the
     reaction energy and uses colors to differentiate them. July 2019 Version.
     Args:

     Kwargs:

     Return:

     """
    G2_obj = nx.Graph()

    colormap_name = ploting_details_dict['colormap']
    colormap = getattr(plt.cm,colormap_name)
    label_dict = {}
    for i in range(len(G1_list)):
        if len(stoich[i]) == 3:
            xpos = stoich[i][0] / stoich[i][1]
            ypos = stoich[i][0]
        elif len(stoich[i]) == 4:
            xpos = stoich[i][0] + stoich[i][1]/ stoich[i][2]
            ypos = stoich[i][1]
        G2_obj.add_node(i,label=G1_list[i].graph["String"],stoich=stoich[i],xpos=xpos,ypos=ypos)
        label_dict[i] = G1_list[i].graph["String"]
        
    if All_models == True:
        reac_idx = [item for sublist in all_reac_idx for item in sublist]
        reac_e = [item for sublist in all_reac_e for item in sublist]
        reac_type = [item for sublist in all_reac_type for item in sublist]
    else:
        reac_idx = all_reac_idx
        reac_e = all_reac_e
        reac_type = all_reac_type

    for i in range(len(reac_idx)):
        ed = [(reac_idx[i][0],reac_idx[i][-1])]
        G2_obj.add_edge(reac_idx[i][0], reac_idx[i][-1], energy=reac_e[i])
        if len(reac_idx[i]) > 2:
            G2_obj.add_edge(reac_idx[i][0], reac_idx[i][-2], energy=reac_e[i])
            ed.append((reac_idx[i][0],reac_idx[i][-2]))

    
    coloredges = [edge[2]["energy"] for edge in G2_obj.edges(data=True)]

    fig = plt.figure(figsize=ploting_details_dict["figsize"],constrained_layout=True)
    ax = fig.subplot_mosaic(mosaic=[["A","B"]],width_ratios=[0.1,1])

    pos = {nd[0]:(nd[1]["xpos"],nd[1]["ypos"]) for nd in G2_obj.nodes(data=True)}

    norm = matplotlib.colors.Normalize(vmin=min(coloredges), vmax=max(coloredges))
    edgecolor = [colormap(norm(e_edge)) for e_edge in coloredges]
    
    nx.draw_networkx_nodes(G2_obj,pos=pos,ax=ax["B"],node_size=50,node_color=ploting_details_dict["node_color"])
    nx.draw_networkx_edges(G2_obj,pos=pos,edge_color=edgecolor,width=2)
    # label_dict = {nd[0]:nd[1] for nd in G2_obj.nodes(data="label")}
    # nx.draw_networkx_labels(G2_obj,pos=pos,labels=label_dict)
    
    col_bar = ax["B"].scatter(coloredges, coloredges, s=0, cmap=colormap, c=coloredges)
    #plt.colorbar(col_bar, ax=ax["A"], fraction=0.01, label='Reaction Energy', location='left')
    plt.colorbar(col_bar, ax=ax["A"],label='Reaction Energy',location="left",fraction=0.33)
    ax["A"].axis("off")
    x_positions = [item for item in nx.get_node_attributes(G2_obj,"xpos").values()]
    y_positions = [item for item in nx.get_node_attributes(G2_obj,"ypos").values()]
    ax["B"].set_xlim([min(x_positions) - 0.1, max(x_positions) + 0.1])
    ax["B"].set_ylim([min(y_positions) - 1, max(y_positions) + 1])

    ax["B"].set_xlabel(ploting_details_dict['x_axis_lab'],fontsize=15,labelpad=20)
    ax["B"].set_ylabel(ploting_details_dict['y_axis_lab'],fontsize=15,labelpad=20)

    fig.suptitle(ploting_details_dict['plot_title'],fontsize=16,y=0.9)

    return fig,ax

def plot_speciation(conc_arr, labels, pH, c0, plot_list=None,ax=None,err_arr=None, col_dict=None,formula_dict=None,
                    raw_concentrations=False,m_idx=0):
    '''Plotting function for speciation diagrams in IPA and HPA systems, supporting the plotting of percentages of different
    species (X and M for XM bimetallic species), raw concentration plots, and error bands
    Args:
        conc_arr: NumPy array of size Nspc x NpH, containing concentration values. With statistical procedures, is often the mean of the
        Nspc x NpH x Nmodel array resulting from speciation calculation.
        labels: list of strings, labels of the species in the diagram.
        pH: array of floats, pH range to compute speciation.
        c0: float, total concentration of the target metal or heteroatom.
        plot_list: list of strings, labels of the species to be plotted. If None, all labels are plotted.
        ax: Matplotlib.Axis object where the diagram will be drawn. If None, a new figure and axis are instantiated.
        err_arr: array of size Nspc x NpH containing the values used to plot error bands, considering a symmetric +- error. Usually,
        the standard deviation of the Nspc x NpH x Nmodel array is used.
        col_dict: dictionary mapping labels to color specifications. If None, tab10 will be used.
        formula_dict: dictionary mapping labels to formatted formula specifications (e.g. generated with Lab_to_Formula). If None,
        raw labels are used for the legend.
        raw_concentrations: boolean. If True, concentrations (in mol/L) will be plotted, without computing percentages for the target element.
        m_idx: integer. Marks the element for which percentages are computed: in XxMmOn HPA species, 0 -> X, 1 -> M.
    Returns:
        ax: Matplotlib.Axis object with the corresponding plot
    '''
    alpha = 0.10
    linewidth = 2.5
    v_ctt2 = [Lab_to_stoich(lab) for lab in labels]
    colors = plt.cm.tab10.colors
    if not plot_list:
        plot_list = labels
    if len(plot_list) > len(colors):
        colors = colors*(int(len(plot_list)/10)+1)
    if col_dict:
        colors = [col_dict[col] for col in plot_list]
    if ax == None:
        fig, ax = plt.subplots(figsize=(5, 5))

    j = 0
    for i, Lab in enumerate(labels):
        if Lab not in plot_list:
            continue
        if formula_dict:
            Lab = formula_dict[Lab]
        multiplier = v_ctt2[i][m_idx] / c0 * 100
        if raw_concentrations:
            multiplier = 1
        percent_concent = conc_arr[i, :] * multiplier
        if i % 2 == 0:
            style = "-"
        else:
            style = "--"

        ax.plot(pH, percent_concent, style, linewidth=linewidth, label=Lab, color=colors[j])

        if err_arr is not None:
            percent_error = err_arr[i, :] * multiplier
            ax.fill_between(pH, percent_concent + percent_error, percent_concent - percent_error,
                            alpha=alpha, color=colors[j])
        j += 1

    ax.set_ylabel("% M")
    ax.set_xlabel("pH")

    return ax


def plot_cluster_means(SuperArr,groups,speciation_labels,pH,c0,
                       col_dict=None,plot_list=None,target_shape=None,m_idx=0):
    '''Convenience function to plot the average speciation for all clusters generated by the SM clusterization pipeline
    Args:
        SuperArr: 3D NumPy array of size Nspc x NpH x Nmodels containing concentration values.
        groups: list of lists of integers, containing indices of the speciation models pertaining to each cluster
        speciation_labels: list of strings, labels of the species for which speciation has been solved
        pH: list of floats, pH values.
        c0: float or list of floats, initial concentration(s) of the reference species.
        col_dict: dictionary mapping species labels to colors
        plot_list: list of strings, species to be included in the plot
        target_shape: tuple of integers, shape of the plot layout. If None, it is decided automatically
        m_idx: integer. Marks the element for which percentages are computed: in XxMmOn HPA species, 0 -> X, 1 -> M.
    Returns:
        fig,axd: Matplotlib figure and dictionary of axis for the figure
    '''
    if not target_shape:
        n_groups = len(groups)
        m = int(np.ceil(np.sqrt(n_groups)))
        if (m*(m-1)) >= n_groups:
            target_shape = (m,m-1)
        else:
            target_shape = (m,m)
    speciation_labels = list(speciation_labels)
    nm = target_shape[0] * target_shape[1]
    exceed_plots = nm - n_groups
    mosaic = np.arange(0, nm).reshape(target_shape)
    # Add another row, for the axis
    mosaic = np.vstack([mosaic, [nm] * target_shape[1]])
    fig = plt.figure(constrained_layout=True, figsize=(target_shape[1]*3, target_shape[0]*2.5))
    axd = fig.subplot_mosaic(mosaic, gridspec_kw={"height_ratios": [1] * target_shape[0] + [0.25]}, sharey=True)
    if not plot_list:
        plot_list = speciation_labels
    for ii, grp in enumerate(groups):
        MiniArr = SuperArr[:, :, np.array(grp)]
        means = np.mean(MiniArr, axis=2)
        stds = np.std(MiniArr, axis=2)
        plot_speciation(means,speciation_labels,pH,c0,plot_list,ax=axd[ii],err_arr=stds,col_dict=col_dict,m_idx=m_idx)
        axd[ii].set_title("Cluster %d (%d)" % (ii, len(grp)))
    for jj in range(ii+1,nm):
        axd[jj].axis("off")

    handles, labels = axd[0].get_legend_handles_labels()
    labels = [Lab_to_Formula(lab) for lab in labels]
    axd[nm].legend(handles, labels, ncol=4, loc="center",fontsize=8)
    axd[nm].axis("off")
    return fig,axd

def plot_boxplot_plot(x_vector, y_position, box_color, box_height, Label=None, ax=None,
                      remove_outliers=False):
    '''Function to easily plot box-and-whisker plots at a requested y-position,
    using data as generated by get_boxplot_data from stats_module.
    Args:
        x_vector: NumPy array of floats, containing values to be used for the boxplot
        y_position: float, vertical position where the boxplot is placed
        box_color: string, hex code for the color of the boxplot (or other color specification supported by matplotlib)
        box_height: float, height of the box
        Label: string, label to tag the boxplot
        ax: Matplotlib axis to place the boxplots at
        remove_outliers: boolean, if True, remove the outliers and regenerate the box-and-whisker plot
    Returns:
        ax: Matplotlib axis
        boxplot: dictionary containing boxplot parameters: quartiles, median, IQR,
                 whisker limits, no. of points in the box, no. of points in the whiskers and no. of outliers
    '''
    boxplot = get_boxplot_data(x_vector)

    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot()

    if remove_outliers:
        x_valid = x_vector[(x_vector < boxplot["wup"]) & (x_vector > boxplot["wdown"])]
        x_vector = x_valid
        boxplot = get_boxplot_data(x_vector)

    ### points (outliers)
    x_outliers = x_vector[(x_vector > boxplot["wup"]) | (x_vector < boxplot["wdown"])]

    y_points = np.tile(y_position, reps=x_outliers.shape[0])
    ax.scatter(x=x_outliers, y=y_points, c="grey", s=2)

    ### box
    box = Rectangle(xy=(boxplot["q1"], y_position - box_height / 2), width=(boxplot["q3"] - boxplot["q1"]),
                    height=box_height, edgecolor="black", facecolor=box_color)
    ax.add_patch(box)

    ### median
    med_line = Line2D(xdata=[boxplot["median"], boxplot["median"]],
                      ydata=[y_position - box_height / 2, y_position + box_height / 2],
                      color="black", solid_capstyle="butt")
    ax.add_line(med_line)

    ### whiskers
    left_horiz = Line2D(xdata=[boxplot["wdown"], boxplot["q1"]], ydata=[y_position, y_position],
                        color="black")
    left_vert = Line2D(xdata=[boxplot["wdown"], boxplot["wdown"]],
                       ydata=[y_position - box_height / 4, y_position + box_height / 4], color="black")
    right_horiz = Line2D(xdata=[boxplot["q3"], boxplot["wup"]], ydata=[y_position, y_position],
                         color="black")
    right_vert = Line2D(xdata=[boxplot["wup"], boxplot["wup"]],
                        ydata=[y_position - box_height / 4, y_position + box_height / 4], color="black")
    [ax.add_line(line) for line in [left_horiz, left_vert, right_horiz, right_vert]]

    ### Annotations
    if Label:
        ax.annotate(Label, (x_vector.max() + 10, y_position), va='center').draggable()

    return ax, boxplot

def plot_best_mods(lgkf_df,sorted_params_df,Kexp,plot_shape=(2,2)):
    '''Auxiliary function to simplify the generation of best speciation models,
    using the information from regressions in lgKf_scaling.
    Args:
        lgkf_df. Pandas DataFrame containing lgkf values for a set of speciation models.
        sorted_params_df. Pandas DataFrame with slope, intercept, regression coefficient, standard error and rmse for each model
        which has been sorted by a parameter of interest (e.g. by rmse to plot best models)
        Kexp. NumPy array with experimental values for rate constants.
        plot_shape. Tuple of integers, marking the no. of panes in the plot.
    Returns:
        fig,ax. Matplotlib Figure and Axis objects.
    '''
    fig, axes = plt.subplots(plot_shape[0], plot_shape[1], sharex=True, sharey=True)
    ax = axes.flatten()
    fig.set_size_inches(4*plot_shape[1],4*plot_shape[0])
    orange = '#e86e00ff'
    blue = '#202252ff'
    Nplots = plot_shape[0]*plot_shape[1]
    for ii in range(Nplots):
        mod_idx = sorted_params_df.index[ii]
        lgkf = lgkf_df.loc[mod_idx,:]
        pars = sorted_params_df.iloc[ii,:]
        m = pars.loc["m"]
        b = pars.loc["b"]
        r2 = (pars.loc["r"])**2
        rmse = pars.loc["rmse"]
        lgkf_scaled = lgkf * m + b
    ##############################################FIGURA 0,0 ###############################################################
        ax[ii].set_ylabel("Experimental Constants " + "$(pK_{f}^{Exp})$", fontsize=12)
        ax[ii].set_xlabel("DFT Constants " + "$(pK_{f}^{DFT})$", fontsize=12)
        ax[ii].plot(lgkf,lgkf_scaled,'-',color=orange)
        ax[ii].plot(lgkf, Kexp, '.', markersize=11, color=blue)
        ax[ii].text(160, 60, "y=%.2fx+(%.2f) and r=%.2f"%(m,b,r2), fontsize=10)
        ax[ii].text(140, 140, "RMSE = %.2f"%rmse, fontsize=10)
    plt.tight_layout()
    return fig,ax

def get_color_phase_diagram(phase_diagram,speciation_labels,col_dict):
    '''Generates an array of RGBA colors to plot phase diagrams according to a color dictionary
    Args:
        phase_diagram: Nvals x NpH array, containing integer indices for the major species at each (C,pH) or (Ratio,pH)
                       point, as produced by phase_diagram_IPA() or phase_diagram_HPA()

        speciation_labels: list of strings, labels of the species for which speciation has been solved
        col_dict: dictionary mapping species labels to colors
    Returns:
        color_array: Nvals x NpH x 4 array containing final RGBA colors to be directly plotted by plt.imshow()
        legend_elements: list of Line2D elements with the colors and labels of the final legend
    '''
    color_list = np.unique(phase_diagram).astype(int)

    if not col_dict:
        nc = len(speciation_labels)
        sample = np.linspace(0,1,nc)
        colors = plt.cm.YlGnBu_r(sample)
        idx_color_dict = {idx:colors[idx] for idx, label in enumerate(speciation_labels)}
        legend_elements = [Line2D([0],[0],marker="o",color=colors[idx],
                                  label=Lab_to_Formula(speciation_labels[idx]),markerfacecolor=colors[idx],markersize=10)
                           for idx in color_list]
    else:
        idx_color_dict = {idx:col_dict[label] for idx, label in enumerate(speciation_labels)}
        legend_elements = [Line2D([0],[0],marker="o",color=col_dict[speciation_labels[idx]],
                                  label=Lab_to_Formula(speciation_labels[idx]),
                                  markerfacecolor=col_dict[speciation_labels[idx]],markersize=10)
                           for idx in color_list]


    rgba_color_dict = {idx: matplotlib.colors.to_rgba(col) for idx, col in idx_color_dict.items()}

    array_tuple = np.vectorize(rgba_color_dict.get)(phase_diagram)
    d1,d2 = phase_diagram.shape
    color_array = np.zeros((d1,d2, 4))
    for jj, layer in enumerate(array_tuple):
        color_array[:, :, jj] = layer

    return color_array,legend_elements
