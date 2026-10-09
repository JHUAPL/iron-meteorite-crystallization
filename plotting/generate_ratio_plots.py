import os
from types import NoneType

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from model.functions import convert_units, make_mixing_lines_ratio, get_meteorite_coords


def generate_ratio_plots(xdata, data1, data2, xopt, opt1, opt2, xunit, unit1, unit2, group, plot_meteorites=True):

    if not os.path.isdir("../../plots/"):
        os.mkdir("../../plots/")

    if not os.path.isdir("../plots/" + group):
        os.mkdir("../plots/" + group)

    two_liq = True

    if unit1 != unit2:
        print('fixing units')
        converted_units = True
        data2 = convert_units(unit2, unit1, data2)
        unit2 = unit1

    elex_sm = xdata[f'SM {xopt}']
    elex_trapped = xdata[f'Trapped SM {xopt}']
    elex_sm_2liq = xdata[f'S2 {xopt}']
    elex_tm1_2liq = xdata[f'Trapped 1 {xopt}']
    elex_tm2_2liq = xdata[f'Trapped 2 {xopt}']
    elex_bulk_trapped = xdata[f'Bulk Trapped {xopt}']

    ele1_sm = data1[f'SM {opt1}']
    ele1_trapped = data1[f'Trapped SM {opt1}']
    ele1_sm_2liq = data1[f'S2 {opt1}']
    ele1_tm1_2liq = data1[f'Trapped 1 {opt1}']
    ele1_tm2_2liq = data1[f'Trapped 2 {opt1}']
    ele1_bulk_trapped = data1[f'Bulk Trapped {opt1}']

    ele2_sm = data2[f'SM {opt2}']
    ele2_trapped = data2[f'Trapped SM {opt2}']
    ele2_sm_2liq_S1 = data2[f'S1 {opt2}']
    ele2_sm_2liq = data2[f'S2 {opt2}']
    ele2_tm1_2liq = data2[f'Trapped 1 {opt2}']
    ele2_tm2_2liq = data2[f'Trapped 2 {opt2}']
    ele2_bulk_trapped = data2[f'Bulk Trapped {opt2}']

    # convert to np arrays for easier indexing
    ele1_tm1_2liq_array = np.array(ele1_tm1_2liq)
    ele1_tm2_2liq_array = np.array(ele1_tm2_2liq)
    ele2_tm1_2liq_array = np.array(ele2_tm1_2liq)
    ele2_tm2_2liq_array = np.array(ele2_tm2_2liq)
    ele1_bulk_trapped_array = np.array(ele1_bulk_trapped)
    ele2_bulk_trapped_array = np.array(ele2_bulk_trapped)
    ele1_sm_2liq_array = np.array(ele1_sm_2liq)
    ele2_sm_2liq_array = np.array(ele2_sm_2liq)

    # # get indices where the trapped concentration is positive
    valid_trapped = np.where((ele1_tm1_2liq_array>0) & (ele1_tm2_2liq_array>0) & (ele2_tm1_2liq_array>0) & (ele2_tm2_2liq_array>0))

    S_start = data1['LM_S'][0]
    P_start = data1['LM_P'][0]
    Ni_start = data1['LM_Ni'][0]
    E1_start = data1[f'LM {opt1}'][0]
    E2_start = data2[f'LM {opt2}'][0]
    Ex_start = xdata[f'LM {xopt}'][0]
    crystallization = np.array(data1['Percent Crystallization'])

    # # MIXING LINE INDICES
    validIdxs = np.where(np.isnan(ele1_sm))

    # one liquid mixing line indices
    # mixing line 1: very top of curve - 0% crystallized
    valid1_1liq = 1
    # mixing line 2: very bottom of one liquid curve (validIdxs[0][1] - 1 = the index right before 2 liquid begins)
    try:
        valid2_1liq = validIdxs[0][1] - 1
    except:
        two_liq = False
        valid2_1liq = len(ele1_sm)-1

    mixing_dict = make_mixing_lines_ratio(xdata, data1, data2, xopt, opt1, opt2)

    colors = sns.color_palette("colorblind")

    plt.rcParams['font.size'] = 12
    plt.rc('axes', titlesize=18)  # Font size for axis titles
    plt.rc('axes', labelsize=18)  # Font size for axis labels
    plt.rc('xtick', labelsize=14)  # Font size for x-axis tick labels
    plt.rc('ytick', labelsize=14)  # Font size for y-axis tick labels
    plt.rc('legend', fontsize=12)
    linesize = 3

    markerSize = 25
    numCols = 5
    log_scale = True

    fig = plt.figure(figsize=(20,10))
    if valid2_1liq > 1:
        # one liquid SM
        plt.plot(elex_sm, ele1_sm/ele2_sm, color=colors[4], linewidth=linesize,label='Solid Metal (One Liquid Field)')
        # one liquid trapped
        plt.plot(elex_trapped, ele1_trapped/ele2_trapped, color=colors[4],  linewidth=linesize, linestyle='dashed', label='Trapped Metal')
        one_liq_lines = mixing_dict['One Liquid']
        for per_cry, opts in one_liq_lines.items():
            ele1_mix = one_liq_lines[per_cry][f'{opt1}']
            ele2_mix = one_liq_lines[per_cry][f'{opt2}']
            plt.plot(ele1_mix, ele2_mix, color=colors[4], linewidth=linesize, linestyle=":")

    if two_liq:
        # two liquid SM
        plt.plot(elex_sm_2liq, ele1_sm_2liq_array/ele2_sm_2liq_array, color=colors[2], linewidth=linesize, label='Solid Metal (Two Liquid Field)')
        # two liquid bulk trapped
        plt.plot(elex_bulk_trapped, ele1_bulk_trapped_array/ele2_bulk_trapped_array, color=colors[2], linewidth=linesize, linestyle='dashdot', label='Bulk Trapped')
        # bulk mixing
        bulk_liq_lines = mixing_dict['Bulk Liquid']
        for per_cry, opts in bulk_liq_lines.items():
            ele1_mix = bulk_liq_lines[per_cry][f'{opt1}']
            ele2_mix = bulk_liq_lines[per_cry][f'{opt2}']
            plt.plot(ele1_mix, ele2_mix, color=colors[2], linewidth=linesize, linestyle=":")


        opac = 0.3
        # S-Rich trapped
        plt.plot(elex_tm1_2liq, ele1_tm1_2liq_array/ele2_tm1_2liq_array, color=colors[1], linewidth=linesize, linestyle='dashdot', alpha=opac, label='Solid from trapped S-rich Liquid')
        # S-Rich mixing
        s_liq_lines = mixing_dict['S-Rich Liquid']
        for per_cry, opts in s_liq_lines.items():
            ele1_mix = s_liq_lines[per_cry][f'{opt1}']
            ele2_mix = s_liq_lines[per_cry][f'{opt2}']
            plt.plot(ele1_mix, ele2_mix, color=colors[1], linewidth=linesize, linestyle=":", alpha=opac)


        # P-Rich trapped
        plt.plot(elex_tm2_2liq, ele1_tm2_2liq_array/ele2_tm2_2liq_array, color=colors[0], linewidth=linesize, linestyle='dashdot', alpha=opac, label='Solid from trapped P-rich Liquid')
        # P-Rich mixing
        p_liq_lines = mixing_dict['P-Rich Liquid']
        for per_cry, opts in p_liq_lines.items():
            ele1_mix = p_liq_lines[per_cry][f'{opt1}']
            ele2_mix = p_liq_lines[per_cry][f'{opt2}']
            plt.plot(ele1_mix, ele2_mix, color=colors[0], linewidth=linesize, linestyle=":", alpha=opac)

    try:
        if two_liq:
            cleaned_opt1 = ele1_sm[~np.isnan(ele1_sm) & ~np.isnan(ele2_sm)].values
            cleaned_opt2 = ele2_sm[~np.isnan(ele1_sm) & ~np.isnan(ele2_sm)].values
            slope = (cleaned_opt2[-1] - cleaned_opt2[0]) / (cleaned_opt1[-1] - cleaned_opt1[0])
            cleaned_opt1_tr = ele1_sm_2liq_array[~np.isnan(ele1_sm_2liq_array) & ~np.isnan(ele2_sm_2liq_array)]
            cleaned_opt2_tr = ele2_sm_2liq_array[~np.isnan(ele1_sm_2liq_array) & ~np.isnan(ele2_sm_2liq_array)]
            slope_tr = (cleaned_opt2_tr[-1] - cleaned_opt2_tr[0]) / (cleaned_opt1_tr[-1] - cleaned_opt1_tr[0])
        else:
            slope = (ele2_sm[-1] - ele2_sm[0]) / (ele1_sm[-1] - ele1_sm[0])
            slope_tr = 0
    except:
        slope = 0
        slope_tr = 0

    ax = plt.gca()
    plot_info = (f'Staring Comps:\n {opt2}: {E2_start:.2f} {unit2}\n {opt1}: {E1_start:.2f} {unit1}\n '
                 f'S: {S_start:.2f} wt%\n P: {P_start:.2f} wt%\n Ni: {Ni_start:.2f} wt%')


    if slope >= 0 or slope_tr >= 0:
        ax.text(0.98, 0.05, plot_info,
                transform=ax.transAxes,
                ha='right',
                va='bottom',
                fontsize=10,
                bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.5))
        # for increasing slope, place text in bottom right
        # ax.annotate(plot_info, (1,0), xycoords='axes fraction', textcoords='offset points',
        #             transform=ax.transAxes, horizontalalignment='right', verticalalignment='bottom')
    else:
        ax.text(0.02, 0.05, plot_info,
                transform=ax.transAxes,
                ha='left',
                va='bottom',
                fontsize=10,
                bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.5))


    if plot_meteorites:
        try:
            optx_names, meteorites_optx, _, ungrpx_names, ungrpx, _ = get_meteorite_coords(xopt, xunit, group)
            opt1_names, meteorites_opt1, _, ungrp1_names, ungrp1, _ = get_meteorite_coords(opt1, unit1, group)
            opt2_names, meteorites_opt2, _, ungrp2_names, ungrp2, ungrp_ref = get_meteorite_coords(opt2, unit2, group)
            plt.scatter(meteorites_optx, meteorites_opt1/meteorites_opt2, color=colors[7], alpha=0.7, marker='^', s=markerSize, label=group)
            if type(ungrp1) != NoneType and type(ungrp2) != NoneType:
                if np.nansum(ungrp1)>0 and np.nansum(ungrp2)>0:
                    plt.scatter(ungrpx, ungrp1/ungrp2, color=colors[9], alpha=0.7, marker='o', s=markerSize*1.2, label='Ungrouped')

        except:
            print(f'{opt1} or {opt2} not available for meteorite group, {group}')
            log_scale = False

    plt.xlabel(f"{xopt} ({xunit})")
    plt.ylabel(f"{opt1}/{opt2}")
    if log_scale:
        plt.xscale('log')
        plt.yscale('log')

    handles, labels = ax.get_legend_handles_labels()

    # ax.legend(handles=handles, labels=labels)
    plt.legend(handles=handles, labels=labels,loc='upper center', bbox_to_anchor=(0.5, 1.12),
           fancybox=True, shadow=True, ncol=numCols)
    # plt.show()

    plotName = f'Ratio_{opt1}{E1_start:.2f}d{opt2}{E2_start:.2f}_{xopt}{Ex_start:.2f}_S{S_start:.2f}_P{P_start:.2f}_Ni{Ni_start:.2f}.png'
    plt.savefig("../plots/" + group + "/" + plotName)
    plt.close()



