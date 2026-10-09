import os

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import numpy as np


def verification_plot(opt, unit, state_data, data, liq_imm_case):

    os.makedirs("../plots/verification/crystallization/", exist_ok=True)
    os.makedirs("../plots/verification/boundary/", exist_ok=True)
    os.makedirs("../plots/verification/trapped/", exist_ok=True)
    os.makedirs("../plots/verification/liquid/", exist_ok=True)

    LMP = state_data['LM_P']
    LMS = state_data['LM_S']
    LMNi = state_data['LM_Ni']
    LMP1 = state_data['LM1_P']
    LMS1 = state_data['LM1_S']
    LMP2 = state_data['LM2_P']
    LMS2 = state_data['LM2_S']
    PC = np.array(state_data['Percent Crystallization'])
    PL = np.array(state_data['Percent Liquid'])
    PS = np.array(state_data['Percent S-Rich'])


    def line_eqn(x1_pts, y1_pts, x2_pts, y2_pts):
        m = (y2_pts - y1_pts)/(x2_pts - x1_pts)
        b = y1_pts - m*x1_pts

        return m, b

    bulkLS = np.array(state_data['Bulk L_S'])
    bulkLP = np.array(state_data['Bulk L_P'])
    pivot_P = np.array(state_data['Pivot_P'])
    pivot_S = np.array(state_data['Pivot_S'])

    two_liq_indices = np.where(~np.isnan(pivot_P))[0]
    pivot_P_twoliq = pivot_P[two_liq_indices]
    pivot_S_twoliq = pivot_S[two_liq_indices]
    bulkLS_twoliq = bulkLS[two_liq_indices]
    bulkLP_twoliq = bulkLP[two_liq_indices]
    # LMP1_twoliq = LMP1[two_liq_indices].values
    # LMP2_twoliq = LMP2[two_liq_indices].values
    # LMS1_twoliq = LMS1[two_liq_indices].values
    # LMS2_twoliq = LMS2[two_liq_indices].values

    LM = data[f'LM {opt}']
    LM1 = data[f'L1 {opt}']
    LM2 = data[f'L2 {opt}']

    TM = data[f'Trapped SM {opt}']
    TM1 = data[f'Trapped 1 {opt}']
    TM2 = data[f'Trapped 2 {opt}']


    S0 = LMS[0]
    P0 = LMP[0]
    Ni0 = LMNi[0]

    comp_string = f'S{S0:.2f}_P{P0:.2f}_Ni{Ni0:.2f}'
    # colors = sns.color_palette('husl', 5)
    # paired_colors = sns.color_palette("Paired")
    colors = sns.color_palette("colorblind")
    uniform_colors = sns.color_palette("rocket")

    plt.rcParams['font.size'] = 12
    markerSize = 20
    numCols = 3
    linesize = 3

    plt.figure(figsize=(12, 12))
    plt.scatter(LMP, LMS, color=colors[3], label='One Liquid')
    plt.scatter(LMP1, LMS1, color=colors[4], label='S-Rich Liquid')
    plt.scatter(LMP2, LMS2, color=colors[5], label='P-Rich Liquid')
    plt.xlabel('wt% P')
    plt.ylabel('wt% S')
    plt.xlim(0, 10)
    plt.ylim(0, 30)
    # plt.tight_layout()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
               fancybox=True, shadow=True, ncol=numCols)
    # plt.show()
    plt.savefig(f"../plots/verification/boundary/SvP_{comp_string}_Case{liq_imm_case}.png")
    plt.close()

    slopes, intercepts = line_eqn(bulkLP_twoliq, bulkLS_twoliq, pivot_P_twoliq, pivot_S_twoliq)

    plt.figure(figsize=(12, 12))
    plt.scatter(LMP, LMS, color=colors[3], label='One Liquid')
    plt.scatter(LMP1, LMS1, color=colors[4], label='S-Rich Liquid')
    plt.scatter(LMP2, LMS2, color=colors[5], label='P-Rich Liquid')
    step_size = round(len(slopes)/20)-1
    x_range = np.linspace(-1, 10, 10)
    for sidx in range(0, len(slopes), step_size):
        plt.plot(x_range, slopes[sidx] * x_range + intercepts[sidx], alpha=0.5, linestyle=":")
        plt.scatter(pivot_P_twoliq[sidx], pivot_S_twoliq[sidx], s=5, c='black')
        plt.scatter(bulkLP_twoliq[sidx], bulkLS_twoliq[sidx], s=5, c='black')
        # plt.scatter(LMP1_twoliq[sidx], LMS1_twoliq[sidx], s=50, c='teal', marker="*")
        # plt.scatter(LMP2_twoliq[sidx], LMS2_twoliq[sidx], s=50, c='teal', marker="*")
    plt.xlabel('wt% P')
    plt.ylabel('wt% S')
    plt.ylim(0, 33)
    # plt.tight_layout()
    # plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
    #            fancybox=True, shadow=True, ncol=numCols)
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    tie_line_entry = mlines.Line2D([], [], color='black', linestyle=':', label='Tie Lines')
    handles.append(tie_line_entry)
    labels.append('Tie Lines')

    # ax.legend(handles=handles, labels=labels)
    plt.legend(handles=handles, labels=labels,loc='upper center', bbox_to_anchor=(0.5, 1.12),
           fancybox=True, shadow=True, ncol=4)
    # plt.show()
    plt.savefig(f"../plots/verification/boundary/SvP_{comp_string}_tie_lines_Case{liq_imm_case}.png")
    plt.close()

    plt.figure(figsize=(12, 12))
    plt.scatter(PC, LMP, color=colors[3], label='One Liquid')
    plt.scatter(PC, LMP1, color=colors[4], label='S-Rich Liquid')
    plt.scatter(PC, LMP2, color=colors[5], label='P-Rich Liquid')
    plt.ylabel('wt% P')
    plt.xlabel('Amount Crystallized')
    # plt.xlim(2, 100)
    # plt.ylim(0.001, 50)
    # plt.tight_layout()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
               fancybox=True, shadow=True, ncol=numCols)
    # plt.show()
    plt.savefig(f"../plots/verification/crystallization/CrystallizationvP_{comp_string}_Case{liq_imm_case}.png")
    plt.close()

    plt.figure(figsize=(12, 12))
    plt.scatter(PC, LMS, color=colors[3], label='One Liquid')
    plt.scatter(PC, LMS1, color=colors[4], label='S-Rich Liquid')
    plt.scatter(PC, LMS2, color=colors[5], label='P-Rich Liquid')
    plt.ylabel('wt% S')
    plt.xlabel('Amount Crystallized')
    # plt.xlim(2, 100)
    # plt.ylim(0.001, 50)
    # plt.tight_layout()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
               fancybox=True, shadow=True, ncol=numCols)
    # plt.show()
    plt.savefig(f"../plots/verification/crystallization/CrystallizationvS_{comp_string}_Case{liq_imm_case}.png")
    plt.close()


    plt.figure(figsize=(12, 12))
    plt.plot(LMS, LM, color=colors[3], linewidth=linesize, label='One Liquid')
    plt.plot(LMS1, LM1, color=colors[4], linewidth=linesize, label='S-Rich Liquid')
    plt.plot(LMS2, LM2, color=colors[5], linewidth=linesize, label='P-Rich Liquid')
    plt.xlabel('wt% S')
    plt.ylabel(f'Liquid Composition of {opt} ({unit})')
    # plt.xlim(2, 100)
    # plt.ylim(0.001, 50)
    # plt.tight_layout()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
               fancybox=True, shadow=True, ncol=numCols)
    # plt.show()
    plt.savefig(f"../plots/verification/liquid/Liq_{opt}vS_{comp_string}_Case{liq_imm_case}.png")
    plt.close()

    plt.figure(figsize=(12, 12))
    plt.plot(LMS, TM, color=colors[3], linewidth=linesize, label='One Liquid TM')
    plt.plot(LMS1, TM1, color=colors[4], linewidth=linesize, label='S-Rich Liquid TM')
    plt.plot(LMS2, TM2, color=colors[5], linewidth=linesize, label='P-Rich Liquid TM')
    plt.xlabel('wt% S')
    plt.ylabel(f'Composition of {opt} ({unit})')
    # plt.xlim(2, 100)
    # plt.ylim(0.001, 50)
    # plt.tight_layout()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
               fancybox=True, shadow=True, ncol=numCols)
    # plt.show()
    plt.savefig(f"../plots/verification/trapped/Trapped_{opt}vS_{comp_string}_Case{liq_imm_case}.png")
    plt.close()

    plt.figure(figsize=(12, 12))
    plt.scatter(PC, PL, color=colors[3], label='% Total Liquid')
    plt.scatter(PC, PL*PS, color=colors[4], label='% S-Rich Liquid')
    plt.scatter(PC, PL*(1-PS), color=colors[5], label='% P-Rich Liquid')
    plt.ylabel('% Liquid')
    plt.xlabel('% Crystallized')
    # plt.xlim(2, 100)
    # plt.ylim(0.001, 50)
    # plt.tight_layout()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
               fancybox=True, shadow=True, ncol=numCols)
    # plt.show()
    plt.savefig(f"../plots/verification/crystallization/CrystallizationvLiq_{comp_string}_Case{liq_imm_case}.png")
    plt.close()

    plt.figure(figsize=(12, 12))
    plt.scatter(PC, PS, color=colors[4], label='% S-Rich Liquid')
    plt.scatter(PC, (1 - PS), color=colors[5], label='% P-Rich Liquid')
    plt.ylabel('% Liquid')
    plt.xlabel('% Crystallized')
    # plt.xlim(2, 100)
    # plt.ylim(0.001, 50)
    # plt.tight_layout()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
               fancybox=True, shadow=True, ncol=numCols)
    # plt.show()
    plt.savefig(f"../plots/verification/crystallization/CrystallizationvLiqSP_{comp_string}_Case{liq_imm_case}.png")
    plt.close()
