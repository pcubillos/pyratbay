import scipy.constants as sc
"""
Constant values used in the pyrat project
"""

# Units processing (conversion to cm):
units = {'A' :1e-8,
         'nm':1e-7,
         'um':1e-4,
         'mm':1e-1,
         'cm':1,
         'm' :1e+2,
         'km':1e+5,
         'au':sc.au*100,
         'pc':sc.parsec*100,
         # Pressure
         'mbar':1e3,
         'bar':1e6,
         # Temperature:
         'kelvin':1
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

# TLI record length:
tlireclen = 26  # Three doubles and one short
dreclen   =  8  # Double byte length
sreclen   =  2  # Short byte length
