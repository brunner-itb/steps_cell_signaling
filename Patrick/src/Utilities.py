from scipy.constants import N_A

def molar_to_molecules(M, Volume):
    """
    Convert molar concentration to number of molecules.

    Args:
        M (float): Molar concentration in mol/L.
        Volume (float): Volume in cubic meters (m^3).

    Returns:
        float: Number of molecules in the given volume at the specified concentration.
    """
    return M * (N_A * Volume * 1e3)  # 1e3 bcs M = mol/l = 1000 mol/m^3