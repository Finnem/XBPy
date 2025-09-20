import math
A, n = 0.986, 16.689
B, m = 1.335, 3.532

# ---- UFF default r* (Å) per element (one sensible default per element) ----
# r* values are taken from UFF (Rappé et al., 1992).
# For elements with multiple UFF types, we choose a common/default:
#   C: C_3 (sp3), N: N_3 (sp3), O: O_3 (sp3), Si: Si3, P: P_3+3, S: S_3+2, etc.
UFF_RSTAR = {
    # Light / organics
    'H': 2.886,        # H_
    'He': 2.362,       # He4+4
    'Li': 2.451,       # Li
    'Be': 2.745,       # Be3+2
    'B': 4.083,        # B_3
    'C': 3.851,        # C_3
    'N': 3.660,        # N_3
    'O': 3.500,        # O_3
    'F': 3.364,        # F_
    'Ne': 3.243,       # Ne4+4

    # Alkali / alkaline earth (neutrals)
    'Na': 2.983,
    'Mg': 3.021,       # Mg3+2
    'Al': 4.499,       # Al3
    'Si': 4.295,       # Si3
    'P': 4.147,        # P_3+3
    'S': 4.035,        # S_3+2
    'Cl': 3.947,
    'Ar': 3.868,       # Ar4+4

    # 4th period selection (use with care for metals)
    'K': 3.812,        # K_
    'Ca': 3.399,       # Ca6+2
    'Sc': 3.295,       # Sc3+3
    'Ti': 3.175,       # Ti3+4
    'V': 3.144,        # V_3+5
    'Cr': 3.023,       # Cr6+3
    'Mn': 2.961,       # Mn6+2
    'Fe': 2.912,       # Fe3+2
    'Co': 2.872,       # Co6+3
    'Ni': 2.834,       # Ni4+2
    'Cu': 3.495,       # Cu3+1
    'Zn': 2.763,       # Zn3+2
    'Ga': 4.383,       # Ga3+3
    'Ge': 4.280,       # Ge3
    'As': 4.230,       # As3+3
    'Se': 4.205,       # Se3+2
    'Br': 4.189,
    'Kr': 4.141,       # Kr4+4

    # 5th period selection
    'Rb': 4.114,
    'Sr': 3.641,       # Sr6+2
    'Y': 3.345,        # Y_3+3
    'Zr': 3.124,       # Zr3+4
    'Nb': 3.165,       # Nb3+5
    'Mo': 3.052,       # Mo6+6
    'Tc': 2.998,       # Tc6+5
    'Ru': 2.963,       # Ru6+2
    'Rh': 2.929,       # Rh6+3
    'Pd': 2.899,       # Pd4+2
    'Ag': 3.148,       # Agl+1
    'Cd': 2.848,       # Cd3+2
    'In': 4.463,       # In3+3
    'Sn': 4.392,       # Sn3
    'Sb': 4.420,       # Sb3+3
    'Te': 4.470,       # Te3+2
    'I': 4.500,
    'Xe': 4.404,       # Xe4+4

    # 6th period picks (few)
    'Cs': 4.517,
    'Ba': 3.703,       # Ba6+2
    'Hf': 3.141,       # Hf3+4
    'Ta': 3.170,       # Ta3+5
    'W': 3.069,        # W_6+6
}

SQRT_2_1_6 = 2.0 ** (1.0/6.0)  # 1.122462048...
def get_sigma_from_uff(symbol: str) -> float:
    """
    Return σ (Å) for a given element symbol using UFF r* via σ = r*/2^{1/6}.
    Raises KeyError if element is not in the table.
    """
    rstar = UFF_RSTAR[symbol]
    return rstar / SQRT_2_1_6

def lj_impact(symbol: str, distance_angstrom: float) -> float:
    """
    Normalized two-power |LJ| impact for a single-element pair (X–X),
    using σ from UFF (Å). Returns a dimensionless value.
    """
    sigma = get_sigma_from_uff(symbol)
    x = sigma / distance_angstrom
    return A * (x ** n) + B * (x ** m)

def lj_impact_pair(sym1: str, sym2: str, distance_angstrom: float) -> float:
    """
    Normalized two-power |LJ| impact for a mixed pair (A–B),
    using Lorentz mixing for σ: σ_AB = (σ_A + σ_B)/2.
    """
    sigma1 = get_sigma_from_uff(sym1)
    sigma2 = get_sigma_from_uff(sym2)
    sigma_ab = 0.5 * (sigma1 + sigma2)
    x = sigma_ab / distance_angstrom
    return A * (x ** n) + B * (x ** m)