import math
import os
import warnings
from typing import Literal

import numpy as np
import pandas as pd
from scipy.differentiate import derivative
from scipy.optimize import fsolve

am_S = 32.065  # g/mol
am_P = 30.974
am_Fe = 55.845
am_Ni = 58.693

# Functions
def crystallization(C_0: float, D: float, f: float):
    '''
    Crystallize the liquid

    Parameters
    ----------
    :param C_0: the weight concentration of the element E in the completely molten core
    :param D: the solid metal-liquid metal weight ratio partition coefficient
    :param f: the fraction of the molten core which solidifies

    Returns
    -------
    :return:
        C_L - the weight concentration of E in the liquid metal following crystallization
        C_S - the weight concentration of E in the crystallized solid metal
    '''

    C_S = (C_0 * D) / (1 - f + f * D)  # from Chabot and Zhang 2021, eq. 1
    C_L = C_S / D  # from Chabot and Zhang 2021, eq. 2

    return C_S, C_L


def trapped_solid(C_L: float, x: float, y: float, D_troilite_E: float, D_schreibersite_E: float):
    '''
    Calculate the concentration of element E in the trapped solid

    Parameters
    ----------
    :param C_L: the weight concentration of the element E in the liquid
    :param x: the fraction of the trapped liquid melt that solidifies to troilite
    :param y: the fraction of the trapped liquid melt that solidifies to schreibersite
    :param D_troilite_E: The solid/liquid partition coefficient for troilite
    :param D_schreibersite_E: The solid/liquid partition coefficient for schreibersite

    Returns
    -------
    :return: C_S_trap - The concentration of E in the trapped solid
    '''


    # C_S_trap = (C_L - x * C_troilite_E - y * C_schreibersite_E )/ (1 - x - y)  # from Chabot and Zhang 2021, eq. 8
    # C_L = x*C_troilite_E + y * C_schreibersite_E + (1-x-y)*C_S_Trap
    # C_S_Trap/C_troilite_E = D1, C_S_Trap/C_schreibersite_E = D2
    # C_L = x*C_S_Trap/D1 + y*C_S_trap/D2 + (1-x-y)*C_S_Trap
    # C_L = C_S_Trap(x/d1 + y/d2 + (1-x-y))
    # C_S_trap = C_L/(x/D_troilite_E + y/D_schreibersite_E  + (1-x-y))
    C_S_trap = C_L /np.nansum([x / D_troilite_E, y / D_schreibersite_E, (1 - x - y)])

    return C_S_trap


def partition_coeff(D_0: float, beta_i: float, FeDomains: float):
    '''
    Calculate the solid metal/liquid metal partition coefficient, D

    Parameters
    ----------
    :param D_0: solid metal-liquid metal partition coefficient in the light element-free Fe-Ni system
    :param beta_i: constant specific to the element E being fit and the light elements, i (S or P in this model)
    :param FeDomains: the fraction of free Fe atoms available in the liquid metal

    Returns
    -------
    :return: D - the partition coefficient
    '''

    with warnings.catch_warnings(record=True) as w:
        D = D_0 * FeDomains ** beta_i  # from Chabot and Zhang 2021, eq. 3
        if len(w) > 0:
            print(w)

    return D


def calc_FeDomains(X_S: float, X_P: float, X_C: float):
    '''
    Calculate the FeDomains, the fraction of free Fe atoms in the liquid

    Parameters
    ----------
    :param X_S: Molar fraction of sulfur
    :param X_P: Molar fraction of phosphorous
    :param X_C: Molar fraction of carbon

    Returns
    -------
    :return: FeDomains - the fraction of free Fe metals atoms available in the liquid metal
    '''

    if X_S != 0 or X_P != 0:
        FeDomain = (1 - 2 * X_S - 4 * X_P) / (1 - X_S - 3 * X_P)  # from Chabot and Zhang 2021, eq. 4
    elif X_C != 0:
        FeDomain = (1 - 4 * X_C) / (1 - 3 * X_C)  # from Chabot et al. 2017, eq. 2
    else:
        FeDomain = 1

    return FeDomain


def calc_beta_sp(X_S: float, X_P: float, beta_S: float, beta_P: float):
    '''
    Calculate the beta sp constant

    Parameters
    ----------
    :param X_S: Molar fraction of sulfur
    :param X_P: Molar fraction of phosphorous
    :param beta_S: constant for sulfur
    :param beta_P: constant for phosphorous

    Returns
    -------
    :return: beta_sp - constant specific to the element E being fit and the light elements, i (S or P in this model)
    '''

    beta_SP = (2 * X_S / (2 * X_S + 4 * X_P)) * beta_S + (
            4 * X_P / (2 * X_S + 4 * X_P)) * beta_P  # from Chabot and Zhang 2021, eq. 5

    return beta_SP


def calc_next_step_size(f_0: float, Lper: float):
    '''
    Calculate the next step size to maintain equal mass steps as the liquid crystallizes

    Parameters
    ----------
    :param f_0: Initial step size
    :param Lper: percent liquid remaining

    Returns
    -------
    :return: f_n - New step size
    '''

    f_n = f_0 / Lper

    return f_n


def calc_percent_liq(L_n_1: float , f_n: float):
    '''
    Calculate the percent liquid remaining

    Parameters
    ----------
    :param L_n_1: Previous LM concentration (L_(n-1))
    :param f_n: Current step size (fraction of molten core that solidifies)

    Returns
    -------
    :return: L_n - Amount of liquid remaining after step
    '''


    L_n = L_n_1 * (1 - f_n)

    return L_n


def calc_X_i(LMS: float, LMP: float, LMNi: float =0):
    '''
    Calculate the molar fraction of S and P in the liquid

    Parameters
    ----------
    LMS: the weight concentration of S in the liquid metal following crystallization
    LMP: the weight concentration of P in the liquid metal following crystallization
    LMNi: the weight concentration of Ni in the liquid metal following crystallization

    Returns
    -------
    X_S, X_P:  Molar fraction of light elements S and P in the liquid metal respectively

    '''


    am_S = 32.065 #g/mol
    am_P = 30.974
    am_Fe = 55.845
    am_Ni = 58.693

    LMFe = 100-LMS-LMP-LMNi

    # amt sulfur (moles) + amt P (moles) + amt Fe (moles) [100 - %S - %P = remaining = %Fe] = total amt liq
    denom = LMS / am_S + LMP / am_P + LMNi/am_Ni + LMFe/ am_Fe
    X_S = LMS / am_S / denom # amt sulfur (moles)/total amt liq (moles) = molar frac of sulfur
    X_P = LMP / am_P / denom

    return X_S, X_P


def phase_boundary(X_P: float, X_Ni: float, X_Fe: float):
    '''
    Get the phase boundary equation at the current S, P, Ni, Fe concentrations in the liquid (from Ulff-Moller (1998))

    Parameters
    ----------
    :param X_P: the molar fraction of P in the liquid
    :param X_Ni: the molar fraction of nickel in the liquid
    :param X_Fe: the molar fraction of iron in the liquid

    Returns
    -------
    :return: phase boundary equation from Ulff-Moller (1998)

    '''
    # equations and valued from Ulff-Moller (1999) https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1945-5100.1998.tb01626.x
    a0 = -4.25
    a1 = -1.39
    a2 = 0.0534
    a3 = 0.147
    a4 = 0.0269
    alpha_p = 1.36
    gamma_0 = 0.04
    gamma_1 = 0.06

    Lx = np.log(4 * alpha_p * X_P + gamma_1 * (X_Ni / (X_Fe + X_Ni)) ** (0.5))
    return np.exp(a0 + a1 * Lx + a2 * Lx ** 2 + a3 * Lx ** 3 + a4 * Lx ** 4 + gamma_0*X_Ni/(X_Ni+X_Fe))


def derive_pivot(denom: float, X_Ni: float, X_Fe: float):
    '''
    Get the pivot point with the current state conditions

    Parameters
    ----------
    :param denom: concentration in the liquid
    :param X_Ni: the molar fraction of nickel in the liquid
    :param X_Fe: the molar fraction of iron in the liquid

    Returns
    -------
    :return: pivot point used to get the tie line intersetcions

    '''


    x = np.linspace(0.001,10,5000)/am_P/denom # X_P

    end1_S = 29/am_S/denom
    y_pb = phase_boundary(x, X_Ni, X_Fe) # X_S
    y_pb_diff = np.abs(y_pb - end1_S)
    end1_P = x[np.where(y_pb_diff == np.nanmin(y_pb_diff))]

    end2_P = 10/am_P/denom
    end2_S = y_pb[-1]

    m = (end1_S - end2_S) / (end1_P - end2_P)
    b = end2_S - m * end2_P

    def static_tie_line(x):
        return m * x + b

    tan_pt_y = 15.2/am_S/denom
    y_tan_diff = np.abs(y_pb - tan_pt_y)
    tan_pt_x = x[np.where(y_tan_diff == np.nanmin(y_tan_diff))]

    df = derivative(phase_boundary, tan_pt_x, args=(X_Ni, X_Fe), step_direction=1).df

    b_tan = tan_pt_y - df * tan_pt_x

    pivot_x = (b_tan - b)/(m - df)
    pivot_y = static_tie_line(pivot_x)

    return pivot_x[0], pivot_y[0]


def calc_Tie_Line_SPNi(bulkLS: float, bulkLP: float, prevL1P: float, prevL2P: float, bulkLNi: float =
0, check_bounds: bool = False):
    '''
    Get the S and P concentration of each liquid as well as the pivot point

    :param bulkLS: The bulk amount of S in the liquid
    :param bulkLP: The bulk amount of P in the liquid
    :param prevL1P: 1 of 2 previous solutions' x-coordinate
    :param prevL2P: 2 of 2 previous solutions' x-coordinate
    :param bulkLNi: The bulk amount of Ni in the liquid
    :param pivot: The tie-line pivot point
    :param check_bounds: False by default to calculate the tie line intersection. True to check if the
    boundary has been crossed into the two liquid field

    :return: if checking bounds, return boolean to indicate one liquid (false) or two liquid (true).
    Otherwise, return boundary intersection points
    '''

    if bulkLNi > 14:
        bulkLNi = 14

    bulkLFe = 100 - bulkLS - bulkLP - bulkLNi
    denom = bulkLS / am_S + bulkLP / am_P + bulkLNi / am_Ni + bulkLFe / am_Fe

    X_Ni = bulkLNi / am_Ni / denom
    X_Fe = bulkLFe / am_Fe / denom
    X_P = bulkLP/am_P/denom
    X_S = bulkLS/am_S/denom

    if check_bounds:

        def F(x,y):
            return y - phase_boundary(x, X_Ni, X_Fe)

        if F(X_P, X_S)<0:
            twoliq = False
        else:
            twoliq=True

        return twoliq

    else:
        pivot_p, pivot_s = derive_pivot(denom, X_Ni, X_Fe)
        # print(f'pivotP, pivotS: ({pivot_p*am_P*denom}, {pivot_s*am_S*denom})')
        m = (X_S - pivot_s) / (X_P - pivot_p)
        b = pivot_s - m * pivot_p

        def y(x):
            return m*x+b

        def F(x):
            return y(x) - phase_boundary(x, X_Ni, X_Fe)

        f = lambda x: F(x)
        i1 = prevL1P/am_P/denom
        i2 = prevL2P/am_P/denom

        with warnings.catch_warnings(record=True) as w:
            points = fsolve(f, [i1, i2])
            if len(w) > 0:
                print(w[0].message)
        neg = False
        if i1 < 0:
            neg = True

        tries = 100
        while math.isclose(points[0], points[1], abs_tol=0.0001) or (points[0] == i1 and points[1] == i2):
            if neg or i1 < 0:
                i1 = 0
            else:
                i1 = i1 - 0.02
            i2 = i2 + 0.02
            points = fsolve(f, [i1, i2])

            tries = tries - 1

            if tries < 0:
                break

        L1P, L2P = sorted(points)  # smaller P = S-rich, larger P = P-rich

        L1S_Tie = y(L1P)
        L2S_Tie = y(L2P)

        return L1P*am_P*denom, L2P*am_P*denom, L1S_Tie*am_S*denom, L2S_Tie*am_S*denom, pivot_p*am_P*denom, pivot_s*am_S*denom


# Calculate trace element partitioning between S-rich and P-rich liquids
def partition_trace_element(Bulk_LM_E: float, D_L1: float, D_L2: float, L1x:float):
    '''
    Partition the trace element into both liquids

    Parameters
    ----------
    :param Bulk_LM_E: Bulk concentration of E in the liquids
    :param D_L1: Partition coefficient of liquid 1
    :param D_L2: Partition coefficient of liquid 2
    :param L1x: Freaction of liquid 1 in the bulk liquid

    Returns
    :return: CL1 - the concentration of E in liquid 1 CL2 - the concentration of E in liquid 2
    -------

    '''
    CL2 = Bulk_LM_E / (L1x * D_L2/D_L1 + 1 - L1x)
    CL1 = CL2*D_L2/D_L1 # eq. 5 in proposal

    return CL1, CL2

lookup_table = pd.read_csv('./data/parameterization_coeff_lookup.csv', index_col='Element')  # from Chabot et al. 2017, Table 2


def coeff_lookup(element: str, var: str):
    '''
    Look up coefficients from literature for a given element

    :param element: Element to look up coefficient for
    :param var: The coefficient to look up
        - D_0
        - beta_S
        - beta_P
        - beta_C
        - D_troilite
        - D_schreibersite

    :return: coeff - the coefficient specified by var for the specified element
    '''

    coeff = lookup_table.loc[element, var]

    if coeff == '-':
        coeff = 0 # []

    return float(coeff)


def get_default_comp(ele: str, group: str, lookup_file: str):
    '''
    Get the default starting concentration from the literature

    Parameters
    ----------
    :param ele: element abbreviation
    :param group: group name
    :param lookup_file: lookup file with starting compositions

    Returns
    -------
    :return: ele_start - the starting concentration of element E for a given group
    ele_unit - the unit for the starting concentration

    '''
    start_comp = pd.read_csv(lookup_file, encoding='cp1252')

    group_row = start_comp.loc[start_comp['Group'] == group]
    ele_start = float(group_row[ele].values[0])

    unit_row = start_comp.loc[start_comp['Group'] == 'Units']
    ele_unit = unit_row[ele].values[0]

    return ele_start, ele_unit


def convert_units(converting_from: str, converting_to: str, data):
    '''
    Convert data to different units

    Parameters
    ----------
    :param converting_from: current units of the data
    :param converting_to: desired units of the data
    :param data: the value(s) to be converted

    Returns
    -------
    :return: the converted value(s)

    '''
    if converting_to != converting_from:
        if converting_to == 'wt%':
            if converting_from == 'ug/g' or converting_from == 'ppm':
                data = data/1e4
            elif converting_from == 'mg/g':
                data = data/10
            elif converting_from == 'ng/g' or  converting_from == 'ppb':
                data = data/1e7
            else:
                print('INVALID UNITS')
        if converting_to == 'mg/g':
            if converting_from == 'ug/g' or converting_from == 'ppm':
                data = data/1e3
            elif converting_from == 'wt%':
                data = data*10
            elif converting_from == 'ng/g' or  converting_from == 'ppb':
                data = data/1e6
            else:
                print('INVALID UNITS')
        elif converting_to == "ppm" or converting_to == 'ug/g':
            if converting_from == 'wt%':
                data = data*1e4
            elif converting_from == 'mg/g':
                data = data * 1e3
            elif converting_from == 'ng/g' or converting_from == 'ppb':
                data = data / 1e3
            else:
                print('INVALID UNITS')

        elif converting_to == 'ng/g' or converting_to == 'ppb':
            if converting_from == 'ug/g' or converting_from == 'ppm':
                data = data * 1e3
            elif converting_from == 'mg/g':
                data = data * 1e6
            elif converting_from == 'wt%':
                data = data * 1e7
            else:
                print('INVALID UNITS')

    return data


def get_meteorite_coords(ele: str, ele_units: str, group:str):
    '''
    Get the meteorite concentrations from literature for a given element and group

    Parameters
    ----------
    :param ele: the element abbreviation
    :param ele_units: the units for the element concentrations
    :param group: the group

    Returns
    -------
    :return: names_grp - the meteorite names
    values_grp - the concentration values for element, ele
    refs_grp - the citation
    names_ungrp - the ungrouped meteorite names thought to be associated with the given group
    values_ungrp - the ungrouped meteorite concentrations
    refs_ungrp - the ungrouped meteorite references
    '''

    meteorites = pd.read_excel('./data/IronMeteorites.xlsx', sheet_name=group, header=0)
    meteorite_names = meteorites['Meteorite'].values[1:].astype(str)
    meteorite_refs = meteorites['Reference'].values[1:].astype(str)
    values = meteorites[ele].values
    meteorite_units = values[0]
    points = values[1:].astype(float)

    conv_points = convert_units(meteorite_units, ele_units, points)

    try:
        ungrpd_index = meteorites.loc[meteorites['Meteorite'] == 'Ungrouped'].index[0]
        values_grp = conv_points[0:ungrpd_index - 1]
        names_grp = meteorite_names[0:ungrpd_index-1]
        refs_grp = meteorite_refs[0:ungrpd_index - 1]
        values_ungrp = conv_points[ungrpd_index:]
        names_ungrp = meteorite_names[ungrpd_index:]
        refs_ungrp = meteorite_refs[ungrpd_index:]
    except:
        values_ungrp = None
        names_ungrp = None
        refs_ungrp = None
        values_grp = conv_points
        names_grp = meteorite_names
        refs_grp = meteorite_refs

    return names_grp, values_grp, refs_grp, names_ungrp,values_ungrp, refs_ungrp


def make_mixing_lines_ratio(xdata, data1, data2, xopt: str, opt1: str, opt2: str, num_lines: int =None, per_crystallized: list=None):
    '''
    Get data to plot mixing lines on ratioed plot

    Parameters
    ----------
    :param xdata: element data for the x axis
    :param data1: element data to be ratioed
    :param data2: element data to be ratioed
    :param xopt: element string for x data
    :param opt1: element string for data 1
    :param opt2: element string for data 2
    :param num_lines: integer number of mixing lines
    :param per_crystallized: list of percents crystallization to plot mixing lines

    Returns
    -------
    :return: dictionary of mixing line data

    '''

    two_liq = True

    ele1_sm = xdata[f'SM {xopt}']
    ele1_trapped = xdata[f'Trapped SM {xopt}']
    ele1_sm_2liq = xdata[f'S2 {xopt}']
    ele1_tm1_2liq = xdata[f'Trapped 1 {xopt}']
    ele1_tm2_2liq = xdata[f'Trapped 2 {xopt}']
    ele1_bulk_trapped = xdata[f'Bulk Trapped {xopt}']

    crystallization = np.array(data1['Percent Crystallization'])

    ele2_sm = np.array(data1[f'SM {opt1}'])/np.array(data2[f'SM {opt2}'])
    ele2_trapped =  np.array(data1[f'Trapped SM {opt1}'])/np.array(data2[f'Trapped SM {opt2}'])
    ele2_sm_2liq = np.array(data1[f'S2 {opt1}'])/np.array(data2[f'S2 {opt2}'])
    ele2_tm1_2liq = np.array(data1[f'Trapped 1 {opt1}'])/np.array(data2[f'Trapped 1 {opt2}'])
    ele2_tm2_2liq = np.array(data1[f'Trapped 2 {opt1}'])/(data2[f'Trapped 2 {opt2}'])
    ele2_bulk_trapped = np.array(data1[f'Bulk Trapped {opt1}'])/(data2[f'Bulk Trapped {opt2}'])

    # convert to np arrays for easier indexing
    ele1_tm1_2liq_array = np.array(ele1_tm1_2liq)
    ele1_tm2_2liq_array = np.array(ele1_tm2_2liq)
    ele2_tm1_2liq_array = np.array(ele2_tm1_2liq)
    ele2_tm2_2liq_array = np.array(ele2_tm2_2liq)

    # get indices where the trapped concentration is positive
    valid_trapped = np.where((ele1_tm1_2liq_array>0) & (ele1_tm2_2liq_array>0) & (ele2_tm1_2liq_array>0) & (ele2_tm2_2liq_array>0))


    # # MIXING LINE INDICES
    validIdxs = np.where(np.isnan(ele1_sm))

    # one liquid mixing line indices
    # mixing line 1: very top of curve - 0% crystallized
    one_liq_top = 1
    # mixing line 2: very bottom of one liquid curve (validIdxs[0][1] - 1 = the index right before 2 liquid begins)
    try:
        one_liq_bottom = validIdxs[0][1] - 1
    except:
        two_liq = False
        one_liq_bottom = len(ele1_sm)-1

    if two_liq:
        # two liquid  mixing line indices
        # mixing line 1: very top of two liquid curve (validIdxs[0][1] + 50 = 50 indices after 2 liquid begins)
        two_liq_top = validIdxs[0][1]
        # mixing line 2: very bottom of the two liquid curve
        try:
            two_liq_bottom = valid_trapped[0][-1]
        except:
            two_liq_bottom = validIdxs[0][-1]
    else:
        two_liq_top = None
        two_liq_bottom = None


    mixing_dict = {}
    mixing_dict['One Liquid'] = {}
    mixing_dict['S-Rich Liquid'] = {}
    mixing_dict['P-Rich Liquid'] = {}
    mixing_dict['Bulk Liquid'] = {}

    mixing_indices = []

    if per_crystallized is not None:
        max_per_crystallized = crystallization[-1]
        for pc in per_crystallized:
            if pc > max_per_crystallized:
                pc = max_per_crystallized
            diff = np.abs(crystallization - pc)
            pc_idx = np.where(diff == np.min(diff))
            mixing_indices.append(int(pc_idx[0][0]))
    elif num_lines is not None:
        ml_indices = np.linspace(0, len(crystallization)-1, num=num_lines)
        ml_int_indices = [int(ml_i) for ml_i in ml_indices]
        mixing_indices.extend(ml_int_indices)
    else:
        mixing_indices.extend([one_liq_top, one_liq_bottom, two_liq_top, two_liq_bottom])


    for mi in mixing_indices:
        if mi is not None:
            per_cry = crystallization[mi]
            if mi <= one_liq_bottom:
                mixing_dict['One Liquid'][f'PC {per_cry:.4f}'] = {}
                mixing_dict['One Liquid'][f'PC {per_cry:.4f}'][f'{opt1}'] = [
                    ele1_sm[mi] * (100 - i) / 100 + ele1_trapped[mi] * i / 100 for i in range(0, 110, 10)]
                mixing_dict['One Liquid'][f'PC {per_cry:.4f}'][f'{opt2}'] = [
                    ele2_sm[mi] * (100 - i) / 100 + ele2_trapped[mi] * i / 100 for i in range(0, 110, 10)]
            else:
                mixing_dict['S-Rich Liquid'][f'PC {per_cry:.4f}'] = {}
                mixing_dict['S-Rich Liquid'][f'PC {per_cry:.4f}'][f'{opt1}'] = [
                    ele1_sm_2liq[mi] * (100 - i) / 100 + ele1_tm1_2liq[mi] * i / 100 for i in range(0, 110, 10)]
                mixing_dict['S-Rich Liquid'][f'PC {per_cry:.4f}'][f'{opt2}'] = [
                    ele2_sm_2liq[mi] * (100 - i) / 100 + ele2_tm1_2liq[mi] * i / 100 for i in range(0, 110, 10)]

                mixing_dict['P-Rich Liquid'][f'PC {per_cry:.4f}'] = {}
                mixing_dict['P-Rich Liquid'][f'PC {per_cry:.4f}'][f'{opt1}'] = [
                    ele1_sm_2liq[mi] * (100 - i) / 100 + ele1_tm2_2liq[mi] * i / 100 for i in range(0, 110, 10)]
                mixing_dict['P-Rich Liquid'][f'PC {per_cry:.4f}'][f'{opt2}'] = [
                    ele2_sm_2liq[mi] * (100 - i) / 100 + ele2_tm2_2liq[mi] * i / 100 for i in range(0, 110, 10)]

                mixing_dict['Bulk Liquid'][f'PC {per_cry:.4f}'] = {}
                mixing_dict['Bulk Liquid'][f'PC {per_cry:.4f}'][f'{opt1}'] = [
                    ele1_sm_2liq[mi] * (100 - i) / 100 + ele1_bulk_trapped[mi] * i / 100 for i in range(0, 110, 10)]
                mixing_dict['Bulk Liquid'][f'PC {per_cry:.4f}'][f'{opt2}'] = [
                    ele2_sm_2liq[mi] * (100 - i) / 100 + ele2_bulk_trapped[mi] * i / 100 for i in range(0, 110, 10)]

    return mixing_dict


def make_mixing_lines(data1, data2, opt1, opt2, num_lines=None, per_crystallized=None):
    '''
    Get data to plot mixing lines on element-element plot

    Parameters
    ----------
    :param data1: element data for x axis
    :param data2: element datafor y axis
    :param opt1: element string for data 1
    :param opt2: element string for data 2
    :param num_lines: integer number of mixing lines
    :param per_crystallized: list of percents crystallization to plot mixing lines
    Returns
    -------
    :return: dictionary of mixing line data

    '''

    two_liq = True

    ele1_sm = data1[f'SM {opt1}']
    ele1_trapped = data1[f'Trapped SM {opt1}']
    ele1_sm_2liq = data1[f'S2 {opt1}']
    ele1_tm1_2liq = data1[f'Trapped 1 {opt1}']
    ele1_tm2_2liq = data1[f'Trapped 2 {opt1}']
    ele1_bulk_trapped = data1[f'Bulk Trapped {opt1}']

    crystallization = np.array(data1['Percent Crystallization'])

    ele2_sm = data2[f'SM {opt2}']
    ele2_trapped = data2[f'Trapped SM {opt2}']
    ele2_sm_2liq = data2[f'S2 {opt2}']
    ele2_tm1_2liq = data2[f'Trapped 1 {opt2}']
    ele2_tm2_2liq = data2[f'Trapped 2 {opt2}']
    ele2_bulk_trapped = data2[f'Bulk Trapped {opt2}']

    # convert to np arrays for easier indexing
    ele1_tm1_2liq_array = np.array(ele1_tm1_2liq)
    ele1_tm2_2liq_array = np.array(ele1_tm2_2liq)
    ele2_tm1_2liq_array = np.array(ele2_tm1_2liq)
    ele2_tm2_2liq_array = np.array(ele2_tm2_2liq)

    # get indices where the trapped concentration is positive
    valid_trapped = np.where((ele1_tm1_2liq_array>0) & (ele1_tm2_2liq_array>0) & (ele2_tm1_2liq_array>0) & (ele2_tm2_2liq_array>0))


    # # MIXING LINE INDICES
    validIdxs = np.where(np.isnan(ele1_sm))

    # one liquid mixing line indices
    # mixing line 1: very top of curve - 0% crystallized
    one_liq_top = 1
    # mixing line 2: very bottom of one liquid curve (validIdxs[0][1] - 1 = the index right before 2 liquid begins)
    try:
        one_liq_bottom = validIdxs[0][1] - 1
    except:
        two_liq = False
        one_liq_bottom = len(ele1_sm)-1

    if two_liq:
        # two liquid  mixing line indices
        # mixing line 1: very top of two liquid curve (validIdxs[0][1] + 50 = 50 indices after 2 liquid begins)
        two_liq_top = validIdxs[0][1]
        # mixing line 2: very bottom of the two liquid curve
        try:
            two_liq_bottom = valid_trapped[0][-1]
        except:
            two_liq_bottom = validIdxs[0][-1]
    else:
        two_liq_top = None
        two_liq_bottom = None


    mixing_dict = {}
    mixing_dict['One Liquid'] = {}
    mixing_dict['S-Rich Liquid'] = {}
    mixing_dict['P-Rich Liquid'] = {}
    mixing_dict['Bulk Liquid'] = {}

    mixing_indices = []

    if per_crystallized is not None:
        max_per_crystallized = crystallization[-1]
        for pc in per_crystallized:
            if pc > max_per_crystallized:
                pc = max_per_crystallized
            diff = np.abs(crystallization - pc)
            pc_idx = np.where(diff == np.min(diff))
            mixing_indices.append(int(pc_idx[0][0]))
    elif num_lines is not None:
        ml_indices = np.linspace(0, len(crystallization)-1, num=num_lines)
        ml_int_indices = [int(ml_i) for ml_i in ml_indices]
        mixing_indices.extend(ml_int_indices)

    else:
        mixing_indices.extend([one_liq_top, one_liq_bottom, two_liq_top, two_liq_bottom])


    for mi in mixing_indices:
        if mi is not None:
            per_cry = crystallization[mi]
            if mi <= one_liq_bottom:
                mixing_dict['One Liquid'][f'PC {per_cry:.4f}'] = {}
                mixing_dict['One Liquid'][f'PC {per_cry:.4f}'][f'{opt1}'] = [
                    ele1_sm[mi] * (100 - i) / 100 + ele1_trapped[mi] * i / 100 for i in range(0, 110, 10)]
                mixing_dict['One Liquid'][f'PC {per_cry:.4f}'][f'{opt2}'] = [
                    ele2_sm[mi] * (100 - i) / 100 + ele2_trapped[mi] * i / 100 for i in range(0, 110, 10)]
            else:
                mixing_dict['S-Rich Liquid'][f'PC {per_cry:.4f}'] = {}
                mixing_dict['S-Rich Liquid'][f'PC {per_cry:.4f}'][f'{opt1}'] = [
                    ele1_sm_2liq[mi] * (100 - i) / 100 + ele1_tm1_2liq[mi] * i / 100 for i in range(0, 110, 10)]
                mixing_dict['S-Rich Liquid'][f'PC {per_cry:.4f}'][f'{opt2}'] = [
                    ele2_sm_2liq[mi] * (100 - i) / 100 + ele2_tm1_2liq[mi] * i / 100 for i in range(0, 110, 10)]

                mixing_dict['P-Rich Liquid'][f'PC {per_cry:.4f}'] = {}
                mixing_dict['P-Rich Liquid'][f'PC {per_cry:.4f}'][f'{opt1}'] = [
                    ele1_sm_2liq[mi] * (100 - i) / 100 + ele1_tm2_2liq[mi] * i / 100 for i in range(0, 110, 10)]
                mixing_dict['P-Rich Liquid'][f'PC {per_cry:.4f}'][f'{opt2}'] = [
                    ele2_sm_2liq[mi] * (100 - i) / 100 + ele2_tm2_2liq[mi] * i / 100 for i in range(0, 110, 10)]

                mixing_dict['Bulk Liquid'][f'PC {per_cry:.4f}'] = {}
                mixing_dict['Bulk Liquid'][f'PC {per_cry:.4f}'][f'{opt1}'] = [
                    ele1_sm_2liq[mi] * (100 - i) / 100 + ele1_bulk_trapped[mi] * i / 100 for i in range(0, 110, 10)]
                mixing_dict['Bulk Liquid'][f'PC {per_cry:.4f}'][f'{opt2}'] = [
                    ele2_sm_2liq[mi] * (100 - i) / 100 + ele2_bulk_trapped[mi] * i / 100 for i in range(0, 110, 10)]

    return mixing_dict


def get_model_path(S_start: float , P_start: float, Ni_start: float, f: float, runCase: Literal[0,1,2,3]):
    '''
    Return string name of output model data file path

    Parameters
    ----------
    :param S_start: Starting concentration of sulfur
    :param P_start: Starting concentration of phosphorous
    :param Ni_start: Starting concentration of nickel
    :param f: initial step size, f
    :param runCase:
        0 = no immiscibility
        1 = equilibrium
        2 = S-rich liquid crystallization
        3 = P-rich liquid crystallization.

    Returns
    -------
    :return: string path to output model excel file

    '''
    if not os.path.isdir("../modelData/"):
        os.mkdir("../modelData/")
    return f'../modelData/S{S_start:.2f}_P{P_start:.2f}_Ni{Ni_start:.2f}_Step{f}_Case{runCase}.xlsx'


def read_sheet(model_path: str, ele: str, comp: float, unit: str):
    '''
    Read the data from an excel sheet

    Parameters
    ----------
    :param model_path: path to the model file
    :param ele: element abbreviation
    :param comp: element starting concentration
    :param unit: unit for element starting concentration

    Returns
    :return: data in the excel sheet
    -------

    '''

    data = pd.read_excel(model_path, sheet_name=f'{ele} {comp:.2f} {unit.replace("/", "_")}', index_col=None)
    return data


def append_to_data_dict(data_dict):
    '''
    Append nan to empty lists in the data dict (force lists to be the same length)

    Parameters
    ----------
    :param data_dict: the data dictionary for the simulation

    Returns
    -------
    :return: the appended data dictionary for the simulation

    '''
    for key in data_dict.keys():
        if len(data_dict[key]) < len(data_dict['Percent Crystallization']):
            data_dict[key].append(math.nan)

    return data_dict


def initialize_data_dict(f0: float, LMS0: float, LMP0: float, LMNi0: float):
    '''
    Initialize the simulation state data dictionary

    Parameters
    ----------
    :param f0: the initial crystallization step size
    :param LMS0: the starting S concentration
    :param LMP0: the starting P concentration
    :param LMNi0: the starting Ni concentration

    Returns
    -------
    :return: the simulation data dictionary

    '''

    f_array = [f0] # step size
    crystallization_array = [0] # amount of core crystallized
    liquid_array = [1] # amount of liquid remaining
    LMS_array = [LMS0] # amount of sulfur in the liquid
    LMP_array = [LMP0] # amount of phosphorous in the liquid
    SMP_array = [math.nan] # amount of phosphorous in the solid
    SMNi_array = [math.nan] # amount of nickel in the solid
    LMNi_array = [LMNi0] # amount of nickel in the liquid
    FeDoms_array = [math.nan] # Fe Domains
    LMS_1_array = [math.nan] # amount of sulfur in the sulfur-rich liquid
    LMP_1_array = [math.nan] # amount of phosphorous in the sulfur-rich liquid
    SMP_1_array = [math.nan] # amount of phosphorous in the solid that forms from the sulfur-rich liquid
    LMNi_1_array = [math.nan] # amount of nickel in the sulfur-rich liquid
    SMNi_1_array = [math.nan] # amount of nickel in the solid that forms from the sulfur-rich liquid
    LMS_2_array = [math.nan] # amount of sulfur in the phosphorous-rich liquid
    LMP_2_array = [math.nan] # amount of phosphorous in the phosphorous-rich liquid
    SMP_2_array = [math.nan] # amount of phosphorous in the solid that forms from the phosphorous-rich liquid
    LMNi_2_array = [math.nan] # amount of Nickel in the phosphorous-rich liquid
    SMNi_2_array = [math.nan] # amount of Nickel in the solid that forms from the phosphorous-rich liquid
    L1_FeDoms_array = [math.nan] # Fe Domains for the S-rich liquid
    L2_FeDoms_array = [math.nan] # Fe Domains for the P-rich liquid
    percent_S_rich = [math.nan] # percentage of the liquid that's sulfur
    bulk_L_S_array = [math.nan] # bulk S amount in the liquids
    bulk_L_P_array = [math.nan] # bulk P amount in the liquids
    bulk_L_Ni_array = [math.nan] # bulk Ni amount in the liquids
    pivot_p_array = [math.nan]
    pivot_s_array = [math.nan]

    data_dict = {'f': f_array, 'Percent Crystallization': crystallization_array, 'Percent Liquid': liquid_array,
                 'LM_S': LMS_array, 'LM_P': LMP_array, 'SM_P': SMP_array,'LM_Ni': LMNi_array, 'SM_Ni': SMNi_array,
                 'FeDomains': FeDoms_array, 'Pivot_P': pivot_p_array, 'Pivot_S': pivot_s_array,
                 'LM1_S': LMS_1_array, 'LM1_P': LMP_1_array, 'SM1_P': SMP_1_array,'LM1_Ni': LMNi_1_array, 'SM1_Ni': SMNi_1_array,
                 'LM2_S': LMS_2_array, 'LM2_P': LMP_2_array, 'SM2_P': SMP_2_array, 'LM2_Ni': LMNi_2_array, 'SM2_Ni': SMNi_2_array,
                 'FeDomains L1': L1_FeDoms_array, 'FeDomains L2': L2_FeDoms_array,
                 'Percent S-Rich': percent_S_rich,
                 'Bulk L_S': bulk_L_S_array, 'Bulk L_P': bulk_L_P_array, 'Bulk L_Ni': bulk_L_Ni_array}

    return data_dict


def initialize_element_data_dict(element: str, LM0: float):
    '''
    Initialize the simulation element state data dictionary

    Parameters
    ----------
    :param element: the element abbreviation
    :param LM0: the starting concentration of element

    Returns
    -------
    :return: the simulation element state data dictionary

    '''

    crystallization_array = [0] # amount of core crystallized
    ele_beta_array = [math.nan] # beta factor for the trace element
    D_array = [math.nan] # partition coefficient for the trace element
    ele_LM_array = [LM0]  # amount of trace element in the solid
    ele_SM_array = [math.nan] # amount of trace element in the solid
    trapped_array = [math.nan] # total amount of trace element trapped
    el_beta_L1_array = [math.nan] # beta factor for the trace element in the S-rich liquid
    el_beta_L2_array = [math.nan] # beta factor for the trace element in the P-rich liquid
    el_D_L1_array = [math.nan] # partition coefficient for the trace element in the S-rich liquid
    el_D_L2_array = [math.nan] # partition coefficient for the trace element in the P-rich liquid
    el_solid1_array = [math.nan] # amount of trace element in the solid that forms from the S-rich liquid
    el_trapped1_array = [math.nan] # amount of trace element trapped from S-rich liquid
    el_trapped2_array = [math.nan] # amount of trace element trapped from P-rich liquid
    el_solid2_array = [math.nan] # amount of trace element in the solid that forms from the P-rich liquid
    el_liq1_array = [math.nan] # amount of trace element in the S-rich liquid
    el_liq2_array = [math.nan] # amount of trace element in the P-rich liquid
    bulk_trapped = [math.nan] # bulk amount of trace element trapped
    bulk_L_el_array = [math.nan] # bulk trace element amount in the liquids

    data_dict = {'Percent Crystallization': crystallization_array,
                 f'beta {element}': ele_beta_array, f'D {element}': D_array,
                 f'LM {element}': ele_LM_array, f'SM {element}': ele_SM_array,
                 f'Trapped SM {element}': trapped_array,
                 f'beta L1 {element}': el_beta_L1_array, f'D L1 {element}': el_D_L1_array,
                 f'L1 {element}': el_liq1_array, f'S1 {element}': el_solid1_array,
                 f'Trapped 1 {element}': el_trapped1_array,
                 f'beta L2 {element}': el_beta_L2_array, f'D L2 {element}': el_D_L2_array,
                 f'L2 {element}': el_liq2_array, f'S2 {element}': el_solid2_array,
                 f'Trapped 2 {element}': el_trapped2_array,
                 f'Bulk Trapped {element}': bulk_trapped,
                 f'Bulk L {element}': bulk_L_el_array}

    return data_dict

