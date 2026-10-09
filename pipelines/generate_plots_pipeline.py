import os.path

from plotting.generate_element_plot import generate_element_plot
from plotting.verification_plots import verification_plot
import pandas as pd
import numpy as np
from model.functions import get_default_comp, get_model_path, read_sheet
from model.run_crystallization import run_crystallization


def verification(ele, S0, P0, Ni0, E0, unit, f0=0.0005, liq_imm_case=1):
    model_path = get_model_path(S0, P0, Ni0, f0, liq_imm_case)
    try:
        state_data = pd.read_excel(model_path, sheet_name='State Data')
        data = read_sheet(model_path, ele, E0, unit)
    except:
        state_data, data = run_crystallization(ele, E0, unit, S0, P0, Ni0, f0, True, liq_imm_case)

    verification_plot(ele, unit, state_data, data, liq_imm_case)

def pipeline(xOpt, opt, start1, start2, unit1, unit2, S_start, P_start, Ni_start, f0, group, liq_imm_case=1):
    model_path = get_model_path(S_start, P_start, Ni_start, f0, liq_imm_case)
    try:
        state_data = pd.read_excel(model_path, sheet_name='State Data')
        data1 = read_sheet(model_path, xOpt, start1, unit1)
    except:
        state_data, data1 = run_crystallization(xOpt, start1, unit1, S_start, P_start, Ni_start, f0, True, liq_imm_case)

    try:
        data2 = read_sheet(model_path, opt, start2, unit2)
    except:
        _, data2 = run_crystallization(opt, start2, unit2, S_start, P_start, Ni_start, f0, True, liq_imm_case)

    generate_element_plot(state_data, data1, data2, xOpt, opt, unit1, unit2, group, liq_imm_case)


if __name__ == '__main__':
    print('start')

    # TODO: need to make arguments for running this on command line
    lookup_file = './data/start_comp_lookup_default.csv'
    # lookup_file = '../data/start_comp_lookup_experiment.csv'
    initial_data = pd.read_csv(lookup_file)
    options = initial_data.columns.values[4:]
    groups = initial_data.iloc[:, 0].values[1:]
    f0 = 0.0005

    # # ============ Specific element and starting comps ============
    # xOpt = 'As'
    # opt = 'Ir'
    # group = 'IIAB'
    # S_start = 15.0  # 9.246153846153847
    # P_start = 0.6  # 0.2830769230769231
    # ele1_start, unit1 = get_start_comp(xOpt, S_start, group, lookup_file)
    # Ni_start, unitNi = get_start_comp('Ni', S_start, group, lookup_file)
    # Ni_start = Ni_start/10 # mg/g to wt%
    # ele2_start, unit2 = get_start_comp(opt, S_start, group, lookup_file)
    # pipeline(xOpt, opt, ele1_start, ele2_start, unit1, unit2, S_start, P_start, Ni_start, group, True)
    # # pipeline('As', 'P', ele1_start, P_start, unit1, 'wt%', S_start, P_start, group, True)
    #
    # verification(xOpt, S_start, P_start, Ni_start, ele1_start, unit1, liq_imm_case=1)


    # # # ============ option loop ============
    # xOpt = 'As'
    # group = 'IIAB'
    # lookup_file = '../data/start_comp_lookup_experiment.csv'
    # S_start, s_unit = get_default_comp('S', group, lookup_file)
    # P_start, P_unit = get_default_comp('P', group, lookup_file)
    #
    # ni_start, ni_unit = get_default_comp('Ni', group, lookup_file)
    # ni_start = ni_start/10
    # co_start, co_unit = get_default_comp('Co', group, lookup_file)
    #
    #
    # ele1_start, unit1 = get_default_comp(xOpt, group, lookup_file)
    #
    # pipeline('As', 'P', ele1_start, P_start, unit1, 'wt%', S_start, P_start, ni_start, f0, group)
    # pipeline('Ni', 'Co', ni_start*10, co_start, ni_unit, co_unit, S_start, P_start, ni_start, f0, group, 1)
    # verification(xOpt, S_start, P_start, ni_start, ele1_start, 'ppm', f0, 1)
    #
    # for opt in options:
    #     if opt == xOpt: #or opt != 'Rh':
    #         continue
    #     print(opt)
    #     ele2_start, unit2 = get_default_comp(opt, group, lookup_file)
    #     if np.isnan(ele2_start):
    #         continue
    #     # verification(opt, S_start, P_start, ele2_start, unit2, f0)
    #     pipeline(xOpt, opt, ele1_start, ele2_start, unit1, unit2, S_start, P_start, ni_start, f0, group, 1)



    # # ============ Experimental single plot ==============
    # lookup_file = './data/start_comp_lookup_experiment.csv'
    # xOpt = 'As'
    # opt = 'Ir'
    # group = 'IIIAB'
    # S_start = 9.0  # 9.246153846153847
    # P_start = 0.3  # 0.2830769230769231
    #
    # ele1_start, unit1 = get_default_comp(xOpt, group, lookup_file)
    # ele2_start, unit2 = get_default_comp(opt, group, lookup_file)
    # ni_start, ni_unit = get_default_comp('Ni', group, lookup_file)
    # pipeline(xOpt, opt, ele1_start, ele2_start, unit1, unit2, S_start, P_start, ni_start/10, 0.0005, group, 0)

    # # ============ group + option loop ============
    f0 = 0.0005
    xOpt = 'As'
    lookup_file = './data/start_comp_lookup_experiment.csv'
    # lookup_file = './start_comp_lookup_default.csv'

    for g in range(0, len(groups)):
        group = groups[g]


        if group != "IIIAB" and group != "IIC":
            continue

        print(group)

        S_start, s_unit = get_default_comp('S', group, lookup_file)
        ele1_start, unit1 = get_default_comp(xOpt, group, lookup_file)
        P_start, P_unit = get_default_comp('P', group, lookup_file)
        ni_start, ni_unit = get_default_comp('Ni', group, lookup_file)
        if ni_unit == "mg/g":
            ni_start = ni_start/10
            ni_unit = 'wt%'
        co_start, co_unit = get_default_comp('Co', group, lookup_file)

        pipeline('Ni', 'Co', ni_start, co_start, ni_unit, co_unit, S_start, P_start, ni_start, f0, group, 1)
        pipeline(xOpt, 'P', ele1_start, P_start, unit1, P_unit, S_start, P_start, ni_start, f0, group, 1)

        # verification(xOpt, S_start, P_start, ni_start, ele1_start, unit1, f0, 1)

        for opt in options:
            print(opt)
            if opt == xOpt:
                continue
            # if opt != "Ir":
            #     continue
            ele2_start, unit2 = get_default_comp(opt, group, lookup_file)
            if np.isnan(ele2_start):
                continue
            pipeline(xOpt, opt, ele1_start, ele2_start, unit1, unit2, S_start, P_start, ni_start, f0, group, 1)
