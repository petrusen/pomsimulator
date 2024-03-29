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
from pomsimulator.modules.text_module import Lab_to_stoich
from pomsimulator.modules.stats_module import get_boxplot_data


def Reaction_Map_3D_monometal(G1_list, all_reac_idx, all_reac_e, all_reac_type, stoich, All_models=True,ploting_details_dict=None):
    """Converts a graph into a reaction map. It weights the edges as a function of the
     reaction energy and uses colors to differentiate them. July 2019 Version."""

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
        y = (len(list_Z) - (list_Z.count(1) + list_Z.count(8))) + list_Z.count(8)
        trivial_zaxis.append(stoich[i][0])
        trivial_yaxis.append(stoich[i][2])
        trivial_xaxis.append(0 + (stoich[i][0]/stoich[i][1]))
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

    G2_obj = nx.Graph()

    colormap_name = ploting_details_dict['colormap']
    colormap = getattr(plt.cm,colormap_name)
    label_dict = {}
    for i in range(len(G1_list)):
        xpos = stoich[i][0] / stoich[i][1] 
        ypos = stoich[i][0]
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
                    raw_concentrations=False):
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
        multiplier = v_ctt2[i][0] / c0 * 100
        if raw_concentrations:
            multiplier = 1
        percent_concent = conc_arr[i, :] * multiplier
        if i % 2 == 0:
            style = "-"
        else:
            style = "--"

        ax.plot(pH, percent_concent, style, linewidth=linewidth, label=Lab, color=colors[j])
        if err_arr != None:
            percent_error = err_arr[i, :] * v_ctt2[i][0] / c0
            ax.fill_between(pH, percent_concent + percent_error, percent_concent - percent_error,
                            alpha=alpha, color=colors[j])
        j += 1

    ax.set_ylabel("% M")
    ax.set_xlabel("pH")

    return ax

def plot_cluster_means(SuperArr,ExpArr,groups,target_shape,speciation_labels,pH,c0,col_dict=None,add_bands=False,
                       plot_list=None,exp_pH = None):

    nm = target_shape[0] * target_shape[1]
    mosaic = np.arange(0, nm).reshape(target_shape)
    # Add another row, for the axis
    mosaic = np.vstack([mosaic, [nm] * target_shape[1]])
    fig = plt.figure(constrained_layout=True, figsize=(target_shape[1]*3, target_shape[0]*2.5))
    axd = fig.subplot_mosaic(mosaic, gridspec_kw={"height_ratios": [1] * target_shape[0] + [0.1]}, sharey=True)
    if not plot_list:
        plot_list = speciation_labels
    if exp_pH is None:
        exp_pH = pH
    plot_speciation(ExpArr,speciation_labels,exp_pH,c0,plot_list,axd[0],err_arr=None,col_dict=col_dict)

    axd[0].set_title("Experimental")
    for ii, grp in enumerate(groups):
        MiniArr = SuperArr[:, :, np.array(grp)]
        means = np.mean(MiniArr, axis=2)
        stds = np.std(MiniArr, axis=2)
        plot_speciation(means,speciation_labels,pH,c0,plot_list,ax=axd[ii+1],err_arr=stds,col_dict=col_dict)
        axd[ii + 1].set_title("Cluster %d (%d)" % (ii, len(grp)))

    handles, labels = axd[0].get_legend_handles_labels()
    axd[nm].legend(handles, labels, ncol=4, loc="center")
    axd[nm].axis("off")
    return fig,axd

def plot_boxplot_plot(x_vector, y_position, box_color, box_height, Label=None, ax=None,
                      remove_outliers=False):
    '''Function to easily plot box-and-whisker plots at a requested y-position,
    using data as generated by get_boxplot_data from stats_module'''
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

