from typing import Literal

import model.functions as functions
from model.chemistry import calc_trapped_solid, calculate_params
from model.state import SimulationState


# L1 = S-rich liquid
# L2 = P-rich liquid

def initialize_state(element: str, LM0: float, LMS0: float, LMP0: float, LMNi0: float, f0: float):
    '''
    Initialize the state data dictionaries

    Parameters
    ----------
    :param element: element abbreviation
    :param LM0: starting concentration of the element
    :param LMS0: starting concentration of S
    :param LMP0: starting concentration of P
    :param LMNi0: starting concentration of Ni
    :param f0: initial crystallization step size

    Returns
    -------
    :return: model state

    '''

    state = SimulationState()

    state.L = 1.0
    state.f = f0
    state.f0 = f0
    state.LM = LM0
    state.LMS = LMS0
    state.LMP = LMP0
    state.LMNi = LMNi0

    state.data = functions.initialize_data_dict(f0, LMS0, LMP0, LMNi0)
    state.eleData = functions.initialize_element_data_dict(element, LM0)

    return state

def run_one_liquid_stage(state, params, element, use_liq_imm=True):
    '''
    Run fractional crystallization while in the one-liquid field

    Parameters
    ----------
    :param state: simulation state
    :param params: simulation constants
    :param element: element abbreviation
    :param use_liq_imm: True/False whether to account for the onset of liquid immiscibility

    '''

    while not state.twoliq and state.L > 1e-10:

        D_E, FeDomains, beta_E = calculate_params(state.LMS, state.LMP, state.LMNi,
                                                  params.element.beta_S, params.element.beta_P, params.element.D0)


        D_P, _, beta_P = calculate_params(state.LMS, state.LMP, state.LMNi,
                                        params.phosphorus.beta_S, params.phosphorus.beta_P, params.phosphorus.D0)

        D_Ni, _, beta_Ni = calculate_params(state.LMS, state.LMP, state.LMNi,
                                            params.nickel.beta_S, params.nickel.beta_P, params.nickel.D0)


        # crystallize liquid E
        SM, state.LM = functions.crystallization(state.LM, D_E, state.f)
        # crystallize liquid S
        _, state.LMS = functions.crystallization(state.LMS, params.D_S, state.f)
        # crystallize liquid P
        SMP, state.LMP = functions.crystallization(state.LMP, D_P, state.f)
        # crystallize liquid Ni
        SMNi, state.LMNi = functions.crystallization(state.LMNi, D_Ni, state.f)

        trapped = calc_trapped_solid(state.LM, state.LMS, state.LMP, params.D_troilite, params.D_schreibersite,
            params.C_FeS_S, params.C_Fe3P_P)

        state.L = functions.calc_percent_liq(state.L, state.f)

        # append values to data dict
        state.data['f'].append( state.f)
        state.data['Percent Crystallization'].append(1 - state.L)
        state.eleData['Percent Crystallization'].append(1 - state.L)
        state.data['Percent Liquid'].append(state.L)

        state.data['LM_S'].append(state.LMS)
        state.data['LM_P'].append(state.LMP)
        state.data['SM_P'].append(SMP)
        state.data['LM_Ni'].append(state.LMNi)
        state.data['SM_Ni'].append(SMNi)
        state.data['Bulk L_S'].append(state.LMS)
        state.data['Bulk L_P'].append(state.LMP)
        state.data['Bulk L_Ni'].append(state.LMNi)
        state.data['FeDomains'].append(FeDomains)

        state.eleData[f'Bulk L {element}'].append(state.LM)
        state.eleData[f'beta {element}'].append(beta_E)
        state.eleData[f'D {element}'].append(D_E)
        state.eleData[f'LM {element}'].append(state.LM)
        state.eleData[f'SM {element}'].append(SM)
        state.eleData[f'Trapped SM {element}'].append(trapped)

        functions.append_to_data_dict(state.data)
        functions.append_to_data_dict(state.eleData)

        # check bounds
        if use_liq_imm:
            state.twoliq = functions.calc_Tie_Line_SPNi(state.LMS, state.LMP, None,None,
                    state.LMNi,True)

        # if monotectic comp is reached, break the loop
        if state.LMS >= 29 or state.LMP >= 10:
            break

        state.f = functions.calc_next_step_size(state.f0,  state.L)

def run_two_liquid_stage(state, params, element, runCase: Literal[0,1,2,3]):
    '''
    Run fractional crystallization while in the two liquid field

    Parameters
    ----------
    :param state: simulation state
    :param params: simulation constants
    :param element: element abbreviation
    :param runCase:
        0 = no immiscibility
        1 = equilibrium
        2 = S-rich liquid crystallization
        3 = P-rich liquid crystallization.

    '''

    L1P = state.LMP - 0.3
    L2P = state.LMP + 0.3
    L1Ni = state.LMNi
    L2Ni = state.LMNi

    if state.LMS>0:
        while state.L > 1e-10:
            # get the S and P liquid concentrations by getting the tie line intersection
            L1P, L2P, L1S, L2S, pivot_p, pivot_s = functions.calc_Tie_Line_SPNi(state.LMS, state.LMP, L1P, L2P, state.LMNi)
            # amount of L1 in the bulk
            L1x = (state.LMS - L2S) / (L1S - L2S)

            D_E_1, FeDomains_1, beta_E_1 = calculate_params(L1S, L1P, L1Ni, params.element.beta_S, params.element.beta_P, params.element.D0)
            D_E_2, FeDomains_2, beta_E_2 = calculate_params(L2S, L2P, L2Ni, params.element.beta_S, params.element.beta_P, params.element.D0)

            D_Ni_1, _, beta_Ni_1 = calculate_params(L1S, L1P, L1Ni, params.nickel.beta_S, params.nickel.beta_P, params.nickel.D0)
            D_Ni_2, _, beta_Ni_2  = calculate_params(L2S, L2P, L2Ni, params.nickel.beta_S, params.nickel.beta_P, params.nickel.D0)

            D_P_1, _, beta_P_1 = calculate_params(L1S, L1P, L1Ni, params.phosphorus.beta_S, params.phosphorus.beta_P, params.phosphorus.D0)
            D_P_2, _, beta_P_2 = calculate_params(L2S, L2P, L2Ni, params.phosphorus.beta_S, params.phosphorus.beta_P, params.phosphorus.D0)

            # Get the element concentration of each liquid
            L1E, L2E = functions.partition_trace_element(state.LM, D_E_1, D_E_2, L1x)
            L1Ni, L2Ni = functions.partition_trace_element(state.LMNi, D_Ni_1, D_Ni_2, L1x)


            # crystallize E
            SM1, LM1 = functions.crystallization(L1E, D_E_1, state.f)
            SM2, LM2 = functions.crystallization(L2E, D_E_2, state.f)

            # crystallize Ni
            Ni_SM1, Ni_LM1 = functions.crystallization(L1Ni, D_Ni_1, state.f)
            Ni_SM2, Ni_LM2 = functions.crystallization(L2Ni, D_Ni_2, state.f)

            # crystallize S
            S_SM1, S_LM1 = functions.crystallization(L1S, params.D_S, state.f)
            S_SM2, S_LM2 = functions.crystallization(L2S, params.D_S, state.f)

            # crystallize P
            P_SM1,  P_LM1 = functions.crystallization(L1P, D_P_1, state.f)
            P_SM2,  P_LM2 = functions.crystallization(L2P, D_P_2, state.f)

            TM1 = calc_trapped_solid(LM1, S_LM1, P_LM1, params.D_troilite, params.D_schreibersite,
            params.C_FeS_S, params.C_Fe3P_P)
            TM2 = calc_trapped_solid(LM2, S_LM2, P_LM2, params.D_troilite, params.D_schreibersite,
            params.C_FeS_S, params.C_Fe3P_P)

            if runCase == 1:
                # recombine liquids into bulk liquid
                state.LMS = S_LM1*L1x + S_LM2*(1-L1x)
                state.LMP = P_LM1*L1x + P_LM2*(1-L1x)
                state.LMNi = Ni_LM1*L1x + Ni_LM2*(1-L1x)
                state.LM = LM1*L1x + LM2*(1-L1x)
                Bulk_trapped = TM1*L1x + TM2*(1-L1x)
            elif runCase == 2:
                state.LMS = S_LM1
                state.LMP = P_LM1
                state.LMNi = Ni_LM1
                state.LM = LM1
                Bulk_trapped = TM1
            elif runCase == 3:
                state.LMS = S_LM2
                state.LMP = P_LM2
                state.LMNi = Ni_LM2
                state.LM = LM2
                Bulk_trapped = TM2

            state.L = functions.calc_percent_liq(state.L, state.f)

            # append to data dict
            state.data['f'].append(state.f)
            state.data['Percent Crystallization'].append(1-state.L)
            state.eleData['Percent Crystallization'].append(1 - state.L)
            state.data['Percent Liquid'].append(state.L)
            state.data['Pivot_P'].append(pivot_p)
            state.data['Pivot_S'].append(pivot_s)
            state.data['LM1_S'].append(S_LM1)
            state.data['LM1_P'].append(P_LM1)
            state.data['SM1_P'].append(P_SM1)
            state.data['LM1_Ni'].append(Ni_LM1)
            state.data['SM1_Ni'].append(Ni_SM1)
            state.data['FeDomains L1'].append(FeDomains_1)
            state.data['LM2_S'].append(S_LM2)
            state.data['LM2_P'].append(P_LM2)
            state.data['SM2_P'].append(P_SM2)
            state.data['LM2_Ni'].append(Ni_LM2)
            state.data['SM2_Ni'].append(Ni_SM2)
            state.data['FeDomains L2'].append(FeDomains_2)
            state.data['Percent S-Rich'].append(L1x)
            state.data['Bulk L_S'].append(state.LMS)
            state.data['Bulk L_P'].append(state.LMP)
            state.data['Bulk L_Ni'].append(state.LMNi)

            state.eleData[f'beta L1 {element}'].append(beta_E_1)
            state.eleData[f'D L1 {element}'].append(D_E_1)
            state.eleData[f'L1 {element}'].append(LM1)
            state.eleData[f'S1 {element}'].append(SM1)
            state.eleData[f'beta L2 {element}'].append(beta_E_2)
            state.eleData[f'D L2 {element}'].append(D_E_2)
            state.eleData[f'L2 {element}'].append(LM2)
            state.eleData[f'S2 {element}'].append(SM2)
            state.eleData[f'Trapped 1 {element}'].append(TM1)
            state.eleData[f'Trapped 2 {element}'].append(TM2)
            state.eleData[f'Bulk Trapped {element}'].append(Bulk_trapped)
            state.eleData[f'Bulk L {element}'].append(state.LM)

            functions.append_to_data_dict(state.data)
            functions.append_to_data_dict(state.eleData)

            # if monotectic comp is reached, break the loop
            if S_LM1 >= 29 or P_LM2 >= 10:
                break # solidify the rest of the core and break the loop
            else:
                state.f = functions.calc_next_step_size(state.f0, state.L)