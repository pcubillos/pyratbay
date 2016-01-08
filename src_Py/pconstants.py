import scipy.constants as sc

"""
Constant values used in the pyrat project.
"""

# Software versioning:

# Pyrat-Bay:
PBAY_VER  = 0   # Major version
PBAY_MIN  = 0   # Minor version
PBAY_REV  = 12  # Revision

# Pyrat version:
PYRAT_VER =  1  # Major version
PYRAT_MIN =  1  # Minor version
PYRAT_REV = 22  # Revision

# Lineread version:
LR_VER    = 6  # Major version
LR_MIN    = 3  # Minor version
LR_REV    = 7  # Revision


# Unit conversion to CGS:
MTC  = 1e-4  # Microns to cm     (MTC = um/cm)
NTC  = 1e-7  # Nanometers to cm  (NTC = nm/cm)

# Convert from eV to cm-1 (kayser):
# planck   = 6.62620e-34  # Planck constant [J * s]
# lumiere  = 2.997925e10  # speed of light  [cm / s]
# electron = 1.602192e-19 # elementary charge [Coulomb]
# kayser2eV = planck * lumiere / electron
eV = 8065.49179

# Distance to cm:
A  = 1e-8  # Angstrom
nm = 1e-7  # Nanometer
um = 1e-4  # Microns
mm = 1e-1  # Millimeter
cm = 1.0   # Centimeter
m  = 1e+2  # Meter
km = 1e+5  # Kilometer
au = sc.au*100      # Astronomical unit
pc = sc.parsec*100  # Parsec
# Pressure:
mbar = 1e3  # Millibar
bar  = 1e6  # Bar
# Temperature:
kelvin = 1.0

# Unified atomic mass to gr:
amu = sc.physical_constants["unified atomic mass unit"][0] * 1e3

# Amagat (Loschmidt number) molecules cm-3:
amagat = sc.physical_constants[
                 "Loschmidt constant (273.15 K, 101.325 kPa)"][0] * 1e-6


# Universal/Physical constants:

# Boltzmann constant in ergs/kelvin:
k = sc.k * 1e7
# Speed of light in cm/s:
c = sc.c * 100.0
# Elementary charge in statcoulombs (from Wolfram Alpha):
e = 4.803205e-10

C1 = 4 * sc.epsilon_0 * sc.m_e * sc.c**2 / sc.e**2 * 0.01  # cm-1
C2 = sc.h * (sc.c * 100.0) / sc.k                          # cm / Kelvin

# String lengths:
maxnamelen = 20
strfmt = "|S%d"%maxnamelen

# TLI record lengths:
tlireclen = 26  # Three doubles and one short
dreclen   =  8  # Double  byte length
ireclen   =  4  # Integer byte length
sreclen   =  2  # Short   byte length
