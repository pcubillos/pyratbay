# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    # Constants
    'h',
    'k',
    'c',
    'G',
    'sigma',
    # Conversion factors
    'eV',
    'A',
    'nm',
    'um',
    'mm',
    'cm',
    'm',
    'km',
    'au',
    'pc',
    'rearth',
    'rjup',
    'rsun',
    'barye',
    'mbar',
    'pascal',
    'bar',
    'atm',
    'gram',
    'kg',
    'mearth',
    'mjup',
    'msun',
    'amu',
    'me',
    'kelvin',
    'sec',
    'min',
    'hour',
    'day',
    'amagat',
    'e',
    'percent',
    'ppt',
    'ppm',
    'none',
    # Other constants
    'C1',
    'C2',
    'C3',
    'tlireclen',
    'dreclen',
    'ireclen',
    'sreclen',
    'ROOT',
    # Choices
    'dbases',
    'rmodes',
    'transmission_rt',
    'emission_rt',
    'rt_paths',
    'retflags',
    'tmodels',
    'chemmodels',
    'radmodels',
    'molmodels',
    'amodels',
    'rmodels',
    'cmodels',
    ]

import os
import scipy.constants as sc

"""
Constant values used in the pyrat project.

Notes
-----
  Solar system constants come from:
  http://nssdc.gsfc.nasa.gov/planetary/factsheet/
"""

# Universal constants in CGS units:
h = sc.h * 1e7  # Planck constant in erg s
k = sc.k * 1e7  # Boltzmann constant in erg K-1
c = sc.c * 1e2  # Speed of light in cm s-1
G = sc.G * 1e3  # Graviational constant in dyne cm2 g-2
sigma = sc.sigma * 1e3  # Stefan-Boltzmann constant in erg s-1 cn-2 K-4

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
rearth = 6.3781e8  # Earth equatorial radius (Prsa et al. 2016)
rjup   = 7.1492e9  # Jupiter equatorial radius (Prsa et al. 2016)
rsun   = 6.957e10  # Sun radius (Prsa et al. 2016)

# Pressure to Barye:
barye  = 1.0    # Barye (CGS units)
mbar   = 1e3    # Millibar
pascal = 1e1    # Pascal (MKS units)
bar    = 1e6    # Bar
atm    = 1.01e6 # Atmosphere

# Mass to grams:
gram   = 1.0           # Gram
kg     = 1.0e3         # Kilogram
mearth = 5.9724e27     # Earth mass
mjup   = 1.8982e30     # Jupiter mass
msun   = 1.9885e33     # Sun mass  (Prsa et al. 2016)
# Unified atomic mass:
amu    = sc.physical_constants['unified atomic mass unit'][0] * 1e3
me     = sc.m_e * 1e3  # Electron mass

# Temperature to Kelvin degree:
kelvin = 1.0

# Time to seconds:
sec = 1.0
min = 60.0
hour = 3600.0
day = 86400.0

# Amagat (Loschmidt number) molecules cm-3:
amagat = sc.physical_constants[
    'Loschmidt constant (273.15 K, 101.325 kPa)'][0] * 1e-6

# Elementary charge in statcoulombs (from Wolfram Alpha):
e = 4.803205e-10

# Other non-physical units:
percent = 1.0e-2  # Percentage
ppt     = 1.0e-3  # Parts per thousand
ppm     = 1.0e-6  # Part per million

# No units:
none = 1

# Other combination of constants:
C1 = 4 * sc.epsilon_0 * sc.m_e * sc.c**2 / sc.e**2 * 0.01  # cm-1
C2 = sc.h * (sc.c * 100.0) / sc.k                          # cm / Kelvin
C3 = sc.pi * e**2 / (me * c**2)                            # cm

# TLI record lengths:
tlireclen = 26  # Three doubles and one short
dreclen   =  8  # Double  byte length
ireclen   =  4  # Integer byte length
sreclen   =  2  # Short   byte length

# Paths:
ROOT = os.path.realpath(os.path.dirname(__file__) + '/../..') + '/'


# Available line-transition databases:
dbases = [
    'Hitran',
    'Exomol',
    'Repack',
    'Pands',
    'Tioschwenke',
    'Voplez',
    'Vald',
    ]

# Running modes:
rmodes = [
    'tli',
    'atmosphere',
    'opacity',
    'spectrum',
    'mcmc',
    ]

# Transmission radiative transfer:
transmission_rt = [
    'transit',
    ]

# Emission radiative transfer:
emission_rt = [
    'emission',
    ]

# Radiative-transfer observing geometry:
rt_paths = transmission_rt + emission_rt

# Retrieval flags:
retflags = [
    'temp',
    'rad',
    'mol',
    'ray',
    'cloud',
    'patchy',
    'mass',
]

# Temperature models:
tmodels = [
   'isothermal',
   'tcea',
   'madhu',
]

# Chemistry models:
chemmodels = [
    'uniform',
    'tea',
]

# Radius-profile models:
radmodels = [
    'hydro_m',
    'hydro_g',
]

# Molecular-abundance models:
molmodels = [
    'vert',
    'scale'
]

# Alkali models:
amodels = [
   'sodium_vdw',
   'potassium_vdw',
]

# Rayleigh models:
rmodels = [
   'dalgarno_H',
   'dalgarno_H2',
   'dalgarno_He',
   'lecavelier',
]

# Cloud/haze models:
cmodels = [
    'deck',
    'ccsgray',
]

