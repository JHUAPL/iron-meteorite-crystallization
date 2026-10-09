import model.functions as functions

def calculate_params(LS, LP, LNi, beta_S,beta_P, D0):
    '''
    Calculate the coefficients needed to crystallize.

    Parameters
    ----------
    :param LS: float
        Sulfur concentration in the liquid
    :param LP: float
        Phosphorous concentration in the liquid
    :param LNi: float
        Nickel concentration in the liquid
    :param beta_S: float
        element beta S constant
    :param beta_P: float
        element beta P constant
    :param D0: float
        partition coefficient in the light element free Fe-Ni system

    Returns
    -------
    :return: D - the partition coefficient
    FeDomains - the fraction of Fe atoms available in the liquid
    beta - the beta SP coefficient

    '''

    X_S, X_P = functions.calc_X_i(LS, LP, LNi)
    FeDomains = functions.calc_FeDomains(X_S, X_P, 0)
    beta = functions.calc_beta_sp(X_S, X_P, beta_S, beta_P)
    D = functions.partition_coeff(D0, beta, FeDomains=FeDomains)

    return D, FeDomains, beta


def calc_trapped_solid(LE, LS, LP, D_troilite, D_schreibersite, C_FeS_S, C_Fe3P_P):
    '''
    Calculate the concentration of the solid that formed from trapped melt

    Parameters
    ----------
    :param LE: float
        concentration of element, E, in the liquid
    :param LS: float
        concentration of S in the liquid
    :param LP: float
        concentration of P in the liquid
    :param D_troilite:
        the solid metal/troilite partition coefficient from Chabot, N. L., Cueva, R. H., Beck, A. W., & Ash, R. D. (2020). Experimental partitioning of trace elements into schreibersite with applications to IIG iron meteorites. Meteoritics & planetary science, 55(4), 726-743.
    :param D_schreibersite: float
        the solid metal/schreibersite partition coefficient from Chabot, N. L., Hamill, C. D., Shread, E. E., Ash, R. D., & Corrigan, C. M. (2025). An experimental study of trace element partitioning into troilite during iron meteorite crystallization. Meteoritics & Planetary Science, 60(5), 1048-1062.
    :param C_FeS_S: float
        The concentration of sulfure in troilite
    :param C_Fe3P_P:
        The concentration of phosphorous in schreibersite

    Returns
    -------
    :return: trapped - the concentration of element, E, in the trapped solid
    '''

    x = LS / C_FeS_S
    y = max(0, (LP + 1.56 * x - 1.56) / (C_Fe3P_P - 1.56))

    trapped = functions.trapped_solid(LE, x, y, D_troilite, D_schreibersite)

    return trapped