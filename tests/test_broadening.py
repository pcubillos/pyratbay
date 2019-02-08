import os
import sys

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay.broadening as pb
import pyratbay.constants  as pc


def test_min_widths():
    min_temp = 500.0     # Kelvin
    min_wn   = 1.0/(0.6*pc.um)
    max_mass = 27.02534  # HCN molecular mass (amu)
    dmin, lmin = pb.min_widths(min_temp, min_wn, max_mass)
    np.testing.assert_approx_equal(dmin, 0.02567273953100574, significant=7)
    np.testing.assert_approx_equal(lmin, 0.00256727395310057, significant=7)


def test_max_widths():
    min_temp =  500.0  # Kelvin
    max_temp = 2500.0  # Kelvin
    max_wn   = 1.0/(12.0*pc.um)
    min_mass = 1.007940  # H molecular mass (amu)
    max_rad  = 2.89/2.0 * pc.A
    max_press = 100 * pc.bar
    dmax, lmax = pb.max_widths(min_temp, max_temp, max_wn, min_mass,
                               max_rad, max_press)
    np.testing.assert_approx_equal(dmax, 0.014862623078508634, significant=7)
    np.testing.assert_approx_equal(lmax, 8.009256827016788,    significant=7)

