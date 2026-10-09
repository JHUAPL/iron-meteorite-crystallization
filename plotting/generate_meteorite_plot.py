import os
from types import NoneType

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model.functions import get_meteorite_coords


def generate_meteorite_plot(groups, x_element, y_element, unitx="ppm", unity="ppm", plot_ungrouped = False):

    if not os.path.isdir("../plots/"):
        os.mkdir("../plots/")

    if groups == 'all':
        xls = pd.ExcelFile('./data/IronMeteorites.xlsx')
        groups = xls.sheet_names

        # lookup_file = './start_comp_lookup_default.csv'
        # initial_data = pd.read_csv(lookup_file)
        # groups = initial_data.iloc[:, 0].values[1:]
        group_list = 'all_groups'
    else:
        group_list = "_".join(groups)


    group_symbol_dict = {"IIC": ["s", "blue"], "IID":["o", "blue"], "IIF":["^", "blue"], "IIIF":["v", "blue"],
                         "IVB": ["D", "blue"], "SBT":["<", "blue"], "IC":[">", "red"], "IIAB":["h", "red"],
                         "IIIAB": ["*", "red"], "IIIE":["p", "red"], "IVA":["X", "red"], "Ungrouped": ["P", "darkgrey"],
                         "IAB": ["1", "grey"], "Grouped Ungrouped":["8", "grey"]}

    # uniform_colors = sns.color_palette("husl",  n_colors=len(groups))
    # uniform_colors = sns.color_palette("flare", n_colors=len(groups))

    plt.rcParams['font.size'] = 12
    plt.rc('axes', titlesize=18)  # Font size for axis titles
    plt.rc('axes', labelsize=18)  # Font size for axis labels
    plt.rc('xtick', labelsize=14)  # Font size for x-axis tick labels
    plt.rc('ytick', labelsize=14)  # Font size for y-axis tick labels
    plt.rc('legend', fontsize=12)
    linesize = 3

    markerSize = 50
    numCols = 6
    log_scale = True

    fig = plt.figure(figsize=(20, 10))

    ax = plt.gca()

    for g in range(0,len(groups)):
        group = groups[g]
        print(group)
        if group == "Antarctic":
            continue
        try:
            if '/' in x_element:
                x_parts = x_element.split('/')
                num_ele = x_parts[0]
                den_ele = x_parts[1]

                num_names, meteorites_num_ele, _, ungrp_num_names, ungrp_num_ele, _ = get_meteorite_coords(num_ele, unitx, group)
                den_names, meteorites_den_ele, _, ungrp_den_names, ungrp_den_ele, _ = get_meteorite_coords(den_ele, unitx, group)

                meteorites_opt1 = meteorites_num_ele/meteorites_den_ele
                if type(ungrp_num_ele) != NoneType and type(ungrp_den_ele) != NoneType:
                    ungrp1 = ungrp_num_ele/ungrp_den_ele
                else:
                    ungrp1 = None
            else:
                opt1_names, meteorites_opt1, _, ungrp1_names, ungrp1, _ = get_meteorite_coords(x_element, unitx, group)

            opt2_names, meteorites_opt2, _ , ungrp2_names, ungrp2, _ = get_meteorite_coords(y_element, unity, group)
            if group != 'IAB':
                plt.scatter(meteorites_opt1, meteorites_opt2, color=group_symbol_dict[group][1],
                            marker=group_symbol_dict[group][0], s=markerSize, label=group)
            if type(ungrp1) != NoneType and type(ungrp2) != NoneType and plot_ungrouped and group != "Ungrouped":
                if np.nansum(ungrp1) > 0 and np.nansum(ungrp2) > 0:
                    if group == "IIIAB":
                        # sHL_idx = np.where(ungrp2_names=='sHL')
                        sHH_idx = np.where(ungrp2_names=='sHH')[0][0]

                        sHL1 = ungrp1[0:sHH_idx-1]
                        sHH1 = ungrp1[sHH_idx:]

                        sHL2 = ungrp2[0:sHH_idx- 1]
                        sHH2 = ungrp2[sHH_idx:]

                        plt.scatter(sHL1, sHL2, color=group_symbol_dict["Grouped Ungrouped"][1],
                                    marker=group_symbol_dict[group][0],
                                    s=markerSize, label=f'IAB sHL')

                        plt.scatter(sHH1, sHH2, color='black',
                                    marker=group_symbol_dict[group][0],
                                    s=markerSize, label=f'IAB sHH')
                    elif group == "IIAB":
                        plt.scatter(ungrp1, ungrp2, color=group_symbol_dict["Grouped Ungrouped"][1],
                                    marker=group_symbol_dict[group][0],
                                    s=markerSize, label=f'IIG')
                    else:
                        plt.scatter(ungrp1, ungrp2, color=group_symbol_dict["Grouped Ungrouped"][1], marker=group_symbol_dict[group][0],
                                    s=markerSize, label=f'{group} Ungrouped')

        except:
            print(f'{x_element} or {y_element} not available for meteorite group, {group}')
            log_scale = False


    if y_element == "Ni" and x_element == "Au" and unity == "mg/g":
        xmin, xmax = ax.get_xlim()
        au_range = np.linspace(xmin, xmax)
        ni_au_ratio = 11.180/0.15 # mg/g / ppm
        # plt.axline((0.15, 11.180), slope=ni_au_ratio, label="Ni/Au Ratio", color="black")
        plt.plot(au_range, au_range*ni_au_ratio, label="Ni/Au Ratio", color = 'black')
        ymin, ymax = ax.get_ylim()
        plt.ylim(10, ymax)

    if y_element == "Sb" and x_element == "Ni/Co":
        plt.axvline(11180/514)
        plt.axhline(157)

    plt.xlabel(x_element + " (" + unitx + ")")
    plt.ylabel(y_element + " (" + unity + ")")
    if log_scale:
        plt.xscale('log')
        plt.yscale('log')
    # plt.xlim(10**-2, 100)
    ymin, ymax = ax.get_ylim()
    # plt.ylim(10, ymax)
    # plt.tight_layout()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
               fancybox=True, shadow=True, ncol=numCols)
    # plt.show()

    x_element = x_element.replace('/', 'd')
    plotName = f'{y_element}_{x_element}.png'
    plt.savefig(f"../plots/{group_list}_{plotName}")
    plt.close()


if __name__ == '__main__':
    generate_meteorite_plot('all', 'Ni/Co', 'Sb', unity="ppb", plot_ungrouped=True)

    generate_meteorite_plot('all', 'Au', 'Ni', unity="mg/g", plot_ungrouped=True)
    generate_meteorite_plot('all', 'Au', 'Sb', unity="ppb", plot_ungrouped=True)
    generate_meteorite_plot('all', 'Ni', 'Ge', unitx="mg/g", plot_ungrouped=True)
    generate_meteorite_plot('all', 'Ni', 'Ir', unitx="mg/g", plot_ungrouped=True)

