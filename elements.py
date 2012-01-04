"""
elements.py

Some standard element properties. Import as global constants.

"""

WEIGHT = {
    "H": 1.00794,
    "He": 4.002602,
    "Li": 6.941,
    "Be": 9.012182,
    "B": 10.811,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984032,
    "Ne": 20.1797,
    "Na": 22.98976928,
    "Mg": 24.3050,
    "Al": 26.9815386,
    "Si": 28.0855,
    "P": 30.973762,
    "S": 32.065,
    "Cl": 35.453,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955912,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938045,
    "Fe": 55.845,
    "Co": 58.933195,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.64,
    "As": 74.92160,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90585,
    "Zr": 91.224,
    "Nb": 92.90638,
    "Mo": 95.96,
    "Tc": 98,
    "Ru": 101.07,
    "Rh": 102.90550,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.60,
    "I": 126.90447,
    "Xe": 131.293,
    "Cs": 132.9054519,
    "Ba": 137.327,
    "La": 138.90547,
    "Ce": 140.116,
    "Pr": 140.90765,
    "Nd": 144.242,
    "Pm": 145,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.92535,
    "Dy": 162.500,
    "Ho": 164.93032,
    "Er": 167.259,
    "Tm": 168.93421,
    "Yb": 173.054,
    "Lu": 174.9668,
    "Hf": 178.49,
    "Ta": 180.94788,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.966569,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.98040,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226,
    "Ac": 227,
    "Th": 232.03806,
    "Pa": 231.03588,
    "U": 238.02891,
    "Np": 237,
    "Pu": 244,
    "Am": 243,
    "Cm": 247,
    "Bk": 247,
    "Cf": 251,
    "Es": 252,
    "Fm": 257,
    "Md": 258,
    "No": 259,
    "Lr": 262,
    "Rf": 265,
    "Db": 268,
    "Sg": 271,
    "Bh": 272,
    "Hs": 270,
    "Mt": 276,
    "Ds": 281,
    "Rg": 280,
    "Cn": 285,
    "Uut": 284,
    "Uuq": 289,
    "Uup": 288,
    "Uuh": 293,
    "Uuo": 294,
    "X": 0,
    "Y": 0,
    "Z": 0}

# If this is a nice list we can just .index or [slice] to get atomic numbers
ATOMIC_NUMBER = [
    "ZERO", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uuo"]

UFF = {
    "H": (2.5711, 0.0440),
    "He": (2.1043, 0.0560),
    "Li": (2.1836, 0.0250),
    "Be": (2.4455, 0.0850),
    "B": (3.6375, 0.1800),
    "C": (3.4309, 0.1050),
    "N": (3.2607, 0.0690),
    "O": (3.1181, 0.0600),
    "F": (2.9970, 0.0500),
    "Ne": (2.8892, 0.0420),
    "Na": (2.6576, 0.0300),
    "Mg": (2.6914, 0.1110),
    "Al": (4.0082, 0.5050),
    "Si": (3.8264, 0.4020),
    "P": (3.6946, 0.3050),
    "S": (3.5948, 0.2740),
    "Cl": (3.5164, 0.2270),
    "Ar": (3.4460, 0.1850),
    "K": (3.3961, 0.0350),
    "Ca": (3.0282, 0.2380),
    "Sc": (2.9355, 0.0190),
    "Ti": (2.8286, 0.0170),
    "V": (2.8010, 0.0160),
    "Cr": (2.6932, 0.0150),
    "Mn": (2.6380, 0.0130),
    "Fe": (2.5943, 0.0130),
    "Co": (2.5587, 0.0140),
    "Ni": (2.5248, 0.0150),
    "Cu": (3.1137, 0.0050),
    "Zn": (2.4616, 0.1240),
    "Ga": (3.9048, 0.4150),
    "Ge": (3.8130, 0.3790),
    "As": (3.7685, 0.3090),
    "Se": (3.7462, 0.2910),
    "Br": (3.7320, 0.2510),
    "Kr": (3.6892, 0.2200),
    "Rb": (3.6652, 0.0400),
    "Sr": (3.2438, 0.2350),
    "Y": (2.9801, 0.0720),
    "Zr": (2.7832, 0.0690),
    "Nb": (2.8197, 0.0590),
    "Mo": (2.7190, 0.0560),
    "Tc": (2.6709, 0.0480),
    "Ru": (2.6397, 0.0560),
    "Rh": (2.6094, 0.0530),
    "Pd": (2.5827, 0.0480),
    "Ag": (2.8045, 0.0360),
    "Cd": (2.5373, 0.2280),
    "In": (3.9761, 0.5990),
    "Sn": (3.9128, 0.5670),
    "Sb": (3.9378, 0.4490),
    "Te": (3.9823, 0.3980),
    "I": (4.0090, 0.3390),
    "Xe": (3.9235, 0.3320),
    "Cs": (4.0242, 0.0450),
    "Ba": (3.2990, 0.3640),
    "La": (3.1377, 0.0170),
    "Ce": (3.1680, 0.0130),
    "Pr": (3.2126, 0.0100),
    "Nd": (3.1850, 0.0100),
    "Pm": (3.1600, 0.0090),
    "Sm": (3.1360, 0.0080),
    "Eu": (3.1119, 0.0080),
    "Gd": (3.0005, 0.0090),
    "Tb": (3.0745, 0.0070),
    "Dy": (3.0540, 0.0070),
    "Ho": (3.0371, 0.0070),
    "Er": (3.0210, 0.0070),
    "Tm": (3.0059, 0.0060),
    "Yb": (2.9890, 0.2280),
    "Lu": (3.2429, 0.0410),
    "Hf": (2.7983, 0.0720),
    "Ta": (2.8241, 0.0810),
    "W": (2.7342, 0.0670),
    "Re": (2.6317, 0.0660),
    "Os": (2.7796, 0.0370),
    "Ir": (2.5302, 0.0730),
    "Pt": (2.4535, 0.0800),
    "Au": (2.9337, 0.0390),
    "Hg": (2.4099, 0.3850),
    "Tl": (3.8727, 0.6800),
    "Pb": (3.8282, 0.6630),
    "Bi": (3.8932, 0.5180),
    "Po": (4.1952, 0.3250),
    "At": (4.2318, 0.2840),
    "Rn": (4.2451, 0.2480),
    "Fr": (4.3654, 0.0500),
    "Ra": (3.2758, 0.4040),
    "Ac": (3.0985, 0.0330),
    "Th": (3.0255, 0.0260),
    "Pa": (3.0504, 0.0220),
    "U": (3.0246, 0.0220),
    "Np": (3.0504, 0.0190),
    "Pu": (3.0504, 0.0160),
    "Am": (3.0121, 0.0140),
    "Cm": (2.9631, 0.0130),
    "Bk": (2.9747, 0.0130),
    "Cf": (2.9515, 0.0130),
    "Es": (2.9391, 0.0120),
    "Fm": (2.9275, 0.0120),
    "Md": (2.9168, 0.0110),
    "No": (2.8936, 0.0110),
    "Lr": (2.8829, 0.0110)
}

VASP_PSEUDO_PREF= {
    "Cu": "Cu_new",
    "Li": "Li_sv",
    "Na": "Na_pv",
    "K": "K_sv",
    "Ca": "Ca_pv",
    "Sr": "Sr_sv",
    "Rb": "Rb_sv",
    "Cs": "Cs_sv",
    "Ba": "Ba_sv",
    "Cu": "Cu_new",
    "Ti": "Ti_pv",
    "V": "V_pv",
    "Cr": "Cr_pv",
    "Mn": "Mn_pv",
    "Sc": "Sc_sv",
    "Y": "Y_sv",
    "Zr": "Zr_sv",
    "Nb": "Nb_pv",
    "Mo": "Mo_pv",
    "Tc": "Tc_pv",
    "Hf": "Hf_pv",
    "Ta": "Ta_pv",
    "W": "W_pv",
    "Os": "Os_pv",
    "Ga": "Ga_d",
    "Ge": "Ge_d",
    "In": "In_d",
    "Sn": "Sn_d",
    "Tl": "Tl_d",
    "Pb": "Pb_d",
    "Bi": "Bi_d",
    "Po": "Po_d",
    "At": "At_d"
}
