# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np
from scipy.ndimage import gaussian_filter1d


def inversion(params, p):
     '''
     Calculates PT profile for inversion case based on Equation (2) from
     Madhusudhan & Seager 2009.

     Parameters
     ----------
     p:  1D float ndarray
         Pressure array, needs to be equally spaced in log space from bottom
         to top of the atmosphere. Must be given in bars.
     params: 1D float ndarray
         Temperature model parameters:
                 a1 - float, exponential factor in Layer 1,
                 empirically determined to be within range (0.2, 0.6).
                 a2 - float, exponential factor in Layer 2,
                 empirically determined to be within range (0.04, 0.5)
                 p1 - floa, pressure boundary between Layer 1 and 2 (in bars).
                 p2 - float, pressure in the middle of tLayer 2
                 p3 - float, pressure boundary between Layers 2 and 3 (in bars).
                 T3 - float, temperature in the Layer 3.

     Returns
     -------
     T_smooth:  1D array of floats, Gaussian smoothed temperatures,
                no kinks on Layer boundaries

     Example
     -------
     # array of pressures, equally spaced in log space
     p = np.logspace(-5, 2, 100)

     # random values
     a1 = np.random.uniform(0.2  , 0.6 )
     a2 = np.random.uniform(0.04 , 0.5 )
     p1 = np.random.uniform(0.001, 0.01)
     p2 = np.random.uniform(0.01 , 1   )
     p3 = np.random.uniform(0.5  , 10  )
     T3 = np.random.uniform(1500 , 1700)
     params = a1, a2, p1, p2, p3, T3

     # calculate temperature
     T_smooth = Inversion(params, p)
     '''
     # Unpack params:
     a1, a2, p1, p2, p3, T3 = params

     # Set p0 (top of the atmosphere):
     p0 = np.amin(p)

     # Calculate temperatures at layer boundaries:
     T2 = T3 - (np.log(p3/p2) / a2)**2
     T0 = T2 + (np.log(p1/p2) / -a2)**2 - (np.log(p1/p0) / a1)**2
     T1 = T0 + (np.log(p1/p0) / a1)**2

     # Defining arrays for every part of the PT profile:
     p_l1     = p[(np.where((p >= min(p)) & (p < p1)))]
     p_l2_pos = p[(np.where((p >= p1)  & (p < p2)))]
     p_l2_neg = p[(np.where((p >= p2)  & (p < p3)))]
     p_l3     = p[(np.where((p >= p3)  & (p <= max(p))))]

     # Layer 1 temperatures:
     T_l1 = (np.log(p_l1/p0) / a1)**2 + T0

     # Layer 2 temperatures (inversion part):
     T_l2_pos = (np.log(p_l2_pos/p2) / -a2)**2 + T2

     # Layer 2 temperatures (decreasing part):
     T_l2_neg = (np.log(p_l2_neg/p2) / a2)**2 + T2

     # Layer 3 temperatures:
     T_l3     = np.linspace(T3, T3, len(p_l3))

     # Concatenating all temperature arrays:
     T_conc = np.concatenate((T_l1, T_l2_pos, T_l2_neg, T_l3))

     # Full PT profile:
     PT_Inver = (T_l1, p_l1, T_l2_pos, p_l2_pos, T_l2_neg,
                 p_l2_neg, T_l3, p_l3, T_conc, T0, T1, T2, T3)

     # Smoothed PT profile:
     sigma = 4
     T_smooth = gaussian_filter1d(T_conc, sigma, mode='nearest')

     return T_smooth


def no_inversion(params, p):
     '''
     Calculates PT profile for inversion case based on Equation (2) from
     Madhusudhan & Seager 2009.

     Parameters
     ----------
     p:  1D float ndarray
         Pressure array, needs to be equally spaced in log space from bottom
         to top of the atmosphere. Must be given in bars.
     params: 1D float ndarray
         Temperature model parameters:
                 a1 - float, exponential factor in Layer 1,
                 empirically determined to be within range (0.2, 0.6).
                 a2 - float, exponential factor in Layer 2,
                 empirically determined to be within range (0.04, 0.5)
                 p1 - floa, pressure boundary between Layer 1 and 2 (in bars).
                 p3 - float, pressure boundary between Layers 2 and 3 (in bars).
                 T3 - float, temperature in the Layer 3.

     Returns
     -------
     T_smooth:  1D array of floats, Gaussian smoothed temperatures,
                no kinks on Layer boundaries

     Example
     -------
     # array of pressures, equally spaced in log space
     p = np.logspace(-5, 2, 100)

     # random values
     a1 = np.random.uniform(0.2  , 0.6 )
     a2 = np.random.uniform(0.04 , 0.5 )
     p1 = np.random.uniform(0.001, 0.01)
     p3 = np.random.uniform(0.5  , 10  )
     T3 = np.random.uniform(1500 , 1700)
     params = a1, a2, p1, p3, T3

     # calculate temperature
     T_smooth = NoInversion(params, p)
     '''
     # Unpack params:
     a1, a2, p1, p3, T3 = params

     # Set p0 (top of the atmosphere):
     p0 = np.amin(p)

     # Calculate temperatures at layer boundaries:
     T1 = T3 - (np.log(p3/p1) / a2)**2.0
     T0 = T1 - (np.log(p1/p0) / a1)**2.0

     # Defining arrays for every part of the PT profile:
     p_l1     = p[np.where((p >= p0) & (p < p1))]
     p_l2_neg = p[np.where((p >= p1) & (p < p3))]
     p_l3     = p[np.where((p >= p3) & (p <= np.amax(p)))]

     # Layer 1 temperatures:
     T_l1 = (np.log(p_l1/p0) / a1)**2 + T0

     # Layer 2 temperatures decreasing part:
     T_l2_neg = (np.log(p_l2_neg/p1) / a2)**2 + T1

     # Layer 3 temperatures:
     T_l3 = np.linspace(T3, T3, len(p_l3))

     # Concatenate all temperature arrays:
     T_conc = np.concatenate((T_l1, T_l2_neg, T_l3))

     # Full PT profile:
     PT_NoInver = (T_l1, p_l1, T_l2_neg, p_l2_neg, T_l3,
                   p_l3, T_conc, T0, T1, T3)

     # Smoothed PT profile:
     sigma = 4
     T_smooth = gaussian_filter1d(T_conc, sigma, mode='nearest')

     return T_smooth
