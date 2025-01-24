from scipy.constants import N_A

def molar_to_molecules(M, Volume):
    # takes Volume in m^3 and returns molecules
    return M * (N_A * Volume * 1e3)  # 1e3 bcs M = mol/l = 1000 mol/m^3