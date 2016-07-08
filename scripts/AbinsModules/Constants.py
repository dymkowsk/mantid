import numpy as np

"""
Constants for instruments and ABINS
"""

# ABINS internal constants

floats_id = np.dtype(np.float64).num
floats_type = np.dtype(np.float64)

DFT_group = "PhononAB" # name of the group in the hdf file in which extracted data from DFT phonon calculations are stored

Q_data_group = "Q_data" # name of the group where Q data is stored

DW_data_group = "DW_factors" # name of the group where Debye-Waller factors are stored

MSQ_data_group = "MSQ" # name of the group where mean square displacements are stored

S_data_group = "S" # name of the group where dynamical factor is stored


all_instruments = ["None", "TOSCA"] # supported instruments

all_sample_forms = ["SingleCrystal", "Powder"] # valid forms of samples

# keywords which define data structure of KpointsData
all_keywords_k_data = ["weight", "value", "frequencies", "atomic_displacements"]

# keywords which define data structure of AtomsData
all_keywords_atoms_data = ["symbol", "fract_coord", "atom", "sort", "mass"]

# symbols of all elements
all_symbols = ["Ac", "Ag", "Al", "Am", "Ar",  "As", "At" , "Au" , "B"  , "Ba", "Be", "Bh", "Bi", "Bk", "Br", "C" , "Ca" ,
               "Cd", "Ce", "Cf", "Cl", "Cm",  "Cn", "Co" , "Cr" , "Cs" , "Cu", "Db", "Ds", "Dy", "Er", "Es", "Eu", "F" ,
               "Fe", "Fl", "Fm", "Fr", "Ga",  "Gd", "Ge" , "H"  , "He" , "Hf", "Hg", "Ho", "Hs", "I" , "In", "Ir", "K" ,
               "Kr", "La", "Li", "Lr", "Lu",  "Lv", "Md" , "Mg" , "Mn" , "Mo", "Mt", "N" , "Na", "Nb", "Nd", "Ne", "Ni",
               "No", "Np", "O" , "Os", "P" ,  "Pa", "Pb" , "Pd" , "Pm" , "Po", "Pr", "Pt", "Pu", "Ra", "Rb", "Re", "Rf",
               "Rg", "Rh", "Rn", "Ru", "S" ,  "Sb", "Sc" , "Se" , "Sg" , "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te",
               "Th", "Ti", "Tl", "Tm", "U" , "Uuo", "Uup", "Uus", "Uut", "V" , "W" , "Xe", "Y" , "Yb", "Zn", "Zr",
               ]


# conversion constants taken from  (http://physics.nist.gov/cgi-bin/cuu/Convert?exp=0&num=1&From=k&To=hr&Action=Only+show+factor)
k_2_hartree = 3.16681050000000e-06 # 1 K * k_2_hartree = 1 Hartree
cm1_2_hartree = 4.556335252767e-6 # 1 cm-1 * cm1_2_hartree = 1 Hartree


# https://en.wikipedia.org/wiki/Atomic_mass_unit
m_2_hartree = 1822.88839  # 1 amu * m2_hartree = 1 Hartree


# Instruments constants
TOSCA_constant = 1 / 16.0 # magic number for TOSCA...

