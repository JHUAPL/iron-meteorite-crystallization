from dataclasses import dataclass
import model.functions as functions


@dataclass
class ElementParameters:
    D0: float
    beta_S: float
    beta_P: float


@dataclass
class SimulationParameters:
    element: ElementParameters
    phosphorus: ElementParameters
    nickel: ElementParameters

    D_troilite: float
    D_schreibersite: float

    C_FeS_S: float = 36.5
    C_Fe3P_P: float = 15.6
    D_S: float = 0.01


def load_parameters(element):

    element_params = ElementParameters(
        D0=functions.coeff_lookup(element, "D_0"),
        beta_S=float(functions.coeff_lookup(element, "beta_S")),
        beta_P=float(functions.coeff_lookup(element, "beta_P"))
    )

    p_params = ElementParameters(
        D0=functions.coeff_lookup("P", "D_0"),
        beta_S=float(functions.coeff_lookup("P", "beta_S")),
        beta_P=float(functions.coeff_lookup("P", "beta_P"))
    )

    ni_params = ElementParameters(
        D0=functions.coeff_lookup("Ni", "D_0"),
        beta_S=float(functions.coeff_lookup("Ni", "beta_S")),
        beta_P=float(functions.coeff_lookup("Ni", "beta_P"))
    )

    return SimulationParameters(
        element=element_params,
        phosphorus=p_params,
        nickel=ni_params,
        D_troilite=float(
            functions.coeff_lookup(element, "D_troilite")
        ),
        D_schreibersite=float(
            functions.coeff_lookup(element, "D_schreibersite")
        )
    )