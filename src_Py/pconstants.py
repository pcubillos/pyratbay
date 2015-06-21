import scipy.constants as sc

"""
Constant values used in the pyrat project.
"""

# Software versioning:

# Pyrat-Bay:
PBAY_VER  = 0  # Major version
PBAY_MIN  = 0  # Minor version
PBAY_REV  = 1  # Revision

# Pyrat version:
PYRAT_VER = 1  # Major version
PYRAT_MIN = 0  # Minor version
PYRAT_REV = 0  # Revision

# Lineread version:
LR_VER    = 6  # Major version
LR_MIN    = 1  # Minor version
LR_REV    = 0  # Revision


MTC  = 1e-4  # Microns to cm     (MTC = um/cm)
NTC  = 1e-7  # Nanometers to cm  (NTC = nm/cm)

# Units processing (conversion to cm):
units = {'A' :1e-8,
         'nm':1e-7,
         'um':1e-4,
         'mm':1e-1,
         'cm':1.0,
         'm' :1e+2,
         'km':1e+5,
         'au':sc.au*100,
         'pc':sc.parsec*100,
         # Pressure
         'mbar':1e3,
         'bar':1e6,
         # Temperature:
         'kelvin':1.0
}

# String lengths:
maxnamelen = 20
strfmt = "|S%d"%maxnamelen

# Unified atomic mass in grams:
u = sc.physical_constants["unified atomic mass unit"][0] * 1e3
# Boltzmann constant in ergs/kelvin:
k = sc.k * 1e7
# Speed of light in cm/s:
c = sc.c * 100.0
# Elementary charge in statcoulombs (from Wolfram Alpha):
e = 4.803205e-10 

C1 = 4 * sc.epsilon_0 * sc.m_e * sc.c**2 / sc.e**2 * 0.01  # cm-1
C2 = sc.h * (sc.c * 100.0) / sc.k                          # cm / Kelvin

# Amagat (Loschmidt number) in mol cm-3:
#amagat = 44.6150e-6
# Amagat (Loschmidt number) in cm-3:
amagat = 2.68678e19

# TLI record length:
tlireclen = 26  # Three doubles and one short
dreclen   =  8  # Double byte length
ireclen   =  4  # Integer byte length
sreclen   =  2  # Short byte length
