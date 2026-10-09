import os
import argparse
from typing import Literal

import pandas as pd

import model.functions as functions
from model.fractional_crystallization_stages import run_one_liquid_stage, initialize_state, run_two_liquid_stage
from model.parameters import load_parameters


def run_crystallization(element: str, LM0: float, unit: Literal["ppm", "ug/g", "mg/g", "wt%", "ng/g", "ppb"],
                        LMS0: float, LMP0: float, LMNi0: float, f0: float=0.0005, save_data: bool=True, run_case: Literal[0,1,2,3] = 1):
    """
    Run the fractional crystallization model.

    Parameters
    ----------
    :param element:
        Element symbol of the trace element.
    :param LM0:
        Initial concentration of the trace element.
    :param unit:
        Unit of LM0.
    :param LMS0:
        Initial sulfur concentration in wt%.
    :param LMP0:
        Initial phosphorus concentration in wt%.
    :param LMNi0:
        Initial nickel concentration in wt%.
    :param f0:
        Initial crystallization step size.
    :param save_data:
        Whether to save model results to Excel.
    :param run_case:
        Crystallization scenario:
        0 = no immiscibility
        1 = equilibrium
        2 = S-rich liquid crystallization
        3 = P-rich liquid crystallization.

    Returns
    -------
    :return: df_state_data: pd.DataFrame of state data (non-trace element concentration evolution data)
    df_ele_state_data: pd.DataFrame of trace element concentration evolution data
    """

    params = load_parameters(element)
    state = initialize_state(element, LM0, LMS0, LMP0, LMNi0, f0)

    if run_case == 0:
        run_one_liquid_stage(state, params, element, False)
    else:
        run_one_liquid_stage(state, params, element, True)

    if state.twoliq:
        run_two_liquid_stage(state, params, element, run_case)

    df_state_data = pd.DataFrame(state.data)
    df_ele_state_data = pd.DataFrame(state.eleData)

    if save_data:
        sheet_exists = False
        model_path = functions.get_model_path(LMS0, LMP0, LMNi0, f0, run_case)
        sheet_name = f'{element} {LM0:.2f} {unit.replace("/", "_")}'
        if not os.path.exists(model_path):
            df_state_data.to_excel(model_path, index=False, sheet_name='State Data')
        else:
            excel_file = pd.ExcelFile(model_path)
            if sheet_name in excel_file.sheet_names:
                sheet_exists = True

        if not sheet_exists:
            with pd.ExcelWriter(model_path, mode='a', engine='openpyxl', if_sheet_exists='new') as writer:
                df_ele_state_data.to_excel(writer, index=False, sheet_name=sheet_name)

    return df_state_data, df_ele_state_data


def main():
    parser = argparse.ArgumentParser(description="Run the crystallization model.")

    # Required arguments
    parser.add_argument("element", help="Element symbol of the trace element (e.g., As).")
    parser.add_argument("LM0", type=float, help="Initial concentration of the trace element.")
    parser.add_argument("unit", choices=["ppm", "ug/g", "mg/g", "wt%", "ng/g", "ppb"], metavar="unit", help="Unit of measurement for the trace element.")
    parser.add_argument("LMS0",type=float,help="Initial sulfur concentration in wt%%.")
    parser.add_argument("LMP0", type=float, help="Initial phosphorus concentration in wt%%.")
    parser.add_argument("LMNi0", type=float,  help="Initial nickel concentration in wt%%.")

    # Optional arguments
    parser.add_argument("--f0",type=float,default=0.0005,help="Initial crystallization step size (default: 0.0005).")
    parser.add_argument("--save-data",action=argparse.BooleanOptionalAction,default=True, help="Save model data to Excel (default: True).")
    parser.add_argument("--run-case",type=int,choices=[0, 1, 2, 3],default=1,
        help=(
            "Crystallization run case: "
            "0=no immiscibility, "
            "1=equilibrium, "
            "2=S-rich crystallization, "
            "3=P-rich crystallization"
            "(default: 1)."
        )
    )


    args = parser.parse_args()
    run_crystallization(args.element, args.LM0, args.unit, args.LMS0, args.LMP0, args.LMNi0, args.f0, args.save_data, args.run_case)


if __name__ == '__main__':
    main()
