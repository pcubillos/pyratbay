# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'ROOT',
    'FILTERS',
    # Other constants
    'tlireclen',
    'dreclen',
    'ireclen',
    'sreclen',
    # Choices
    'dbases',
    'rmodes',
    'samplers',
    'statistics',
    'transmission_rt',
    'emission_rt',
    'rt_paths',
    'retflags',
    'tmodels',
    'chemmodels',
    'radmodels',
    'amodels',
    'rmodels',
    'cmodels',
    'h_ion_models',
]

import os

"""
Constant values used in the pyrat project.
"""

# TLI record lengths:
tlireclen = 26  # Three doubles and one short
dreclen = 8  # Double byte length
ireclen = 4  # Integer byte length
sreclen = 2  # Short byte length

# Paths:
ROOT = os.path.realpath(os.path.dirname(__file__) + '/../..') + '/'
FILTERS = f'{ROOT}pyratbay/data/filters/'


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
    'radeq',
    'retrieval',
]

# Bayesian retrieval algorithms
samplers = [
    'snooker',
    'multinest',
]

# Statistics for posterior plots
statistics = [
    'med_central',
    'max_central',
    'max_like',
    'global_max_like',
]

# Transmission radiative transfer:
transmission_rt = [
    'transit',
]

# Emission radiative transfer:
emission_rt = [
    'emission',
    'emission_two_stream',
]

# Radiative-transfer observing geometry:
rt_paths = transmission_rt + emission_rt

# Retrieval flags:
retflags = [
    'temp',
    'rad',
    'press',
    'mol',
    'ray',
    'cloud',
    'patchy',
    'mass',
    'tstar',
]

# Temperature models:
tmodels = [
   'isothermal',
   'guillot',
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

# H- opacity models:
h_ion_models = [
    'h_ion_john1988',
]
