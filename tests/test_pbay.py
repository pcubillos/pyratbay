# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import pytest
import subprocess

import numpy as np

from conftest import make_config

import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.io as io
from pyratbay.constants import ROOT

os.chdir(ROOT + 'tests')


expected_pressure    = np.logspace(0, 8, 81)
expected_temperature = np.array([
    1046.89433798, 1046.89534525, 1046.89661946, 1046.89823135,
    1046.90027034, 1046.90284953, 1046.90611195, 1046.91023851,
    1046.91545794, 1046.9220595 , 1046.93040895, 1046.94096875,
    1046.95432364, 1046.97121287, 1046.99257099, 1047.01957934,
    1047.05373107, 1047.09691321, 1047.15151021, 1047.22053447,
    1047.30779086, 1047.41808387, 1047.55747799, 1047.73362488,
    1047.95617334, 1048.23728226, 1048.59226036, 1049.04036124,
    1049.60576707, 1050.31879834, 1051.21739055, 1052.34887907,
    1053.77212894, 1055.56003335, 1057.80237786, 1060.60902021,
    1064.11326045, 1068.47516244, 1073.88443242, 1080.56226117,
    1088.76131307, 1098.76283701, 1110.86975475, 1125.39465039,
    1142.64193973, 1162.88419029, 1186.3335257 , 1213.11007182,
    1243.21017731, 1276.47740788, 1312.57904127, 1350.99026501,
    1390.98801603, 1431.65689115, 1471.91099843, 1510.53770339,
    1546.27096545, 1577.9014716 , 1604.42506087, 1625.21631993,
    1640.18896334, 1649.87550785, 1655.34982468, 1657.96419233,
    1658.9854222 , 1659.31327079, 1659.42226561, 1659.49173767,
    1659.56937861, 1659.66627037, 1659.78818839, 1659.94163515,
    1660.13475269, 1660.37777747, 1660.68357589, 1661.06831324,
    1661.55228907, 1662.16097781, 1662.92632193, 1663.88833303,
    1665.09706555])

expected_abundance = np.array([
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     8.6137e-15, 4.5952e-04, 2.4086e-07, 2.6291e-12, 1.4807e-14, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     1.3648e-14, 4.5952e-04, 2.4085e-07, 3.3096e-12, 1.8641e-14, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     2.1632e-14, 4.5952e-04, 2.4085e-07, 4.1667e-12, 2.3468e-14, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     3.4285e-14, 4.5952e-04, 2.4085e-07, 5.2456e-12, 2.9545e-14, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     5.4337e-14, 4.5952e-04, 2.4085e-07, 6.6038e-12, 3.7194e-14, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     8.6118e-14, 4.5952e-04, 2.4085e-07, 8.3137e-12, 4.6825e-14, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     1.3645e-13, 4.5952e-04, 2.4084e-07, 1.0466e-11, 5.8948e-14, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     2.1627e-13, 4.5952e-04, 2.4084e-07, 1.3176e-11, 7.4212e-14, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     3.4267e-13, 4.5952e-04, 2.4083e-07, 1.6586e-11, 9.3427e-14, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     5.4310e-13, 4.5952e-04, 2.4083e-07, 2.0880e-11, 1.1762e-13, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     8.6054e-13, 4.5952e-04, 2.4082e-07, 2.6285e-11, 1.4807e-13, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     1.3635e-12, 4.5952e-04, 2.4081e-07, 3.3089e-11, 1.8641e-13, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     2.1605e-12, 4.5952e-04, 2.4080e-07, 4.1655e-11, 2.3468e-13, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     3.4226e-12, 4.5952e-04, 2.4078e-07, 5.2435e-11, 2.9544e-13, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     5.4217e-12, 4.5952e-04, 2.4077e-07, 6.6002e-11, 3.7193e-13, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     8.5864e-12, 4.5952e-04, 2.4074e-07, 8.3077e-11, 4.6823e-13, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     1.3598e-11, 4.5952e-04, 2.4071e-07, 1.0457e-10, 5.8947e-13, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     2.1525e-11, 4.5952e-04, 2.4067e-07, 1.3161e-10, 7.4208e-13, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     3.4073e-11, 4.5952e-04, 2.4062e-07, 1.6563e-10, 9.3421e-13, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     5.3909e-11, 4.5952e-04, 2.4056e-07, 2.0843e-10, 1.1761e-12, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     8.5250e-11, 4.5952e-04, 2.4047e-07, 2.6225e-10, 1.4806e-12, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     1.3474e-10, 4.5952e-04, 2.4037e-07, 3.2993e-10, 1.8638e-12, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     2.1282e-10, 4.5952e-04, 2.4025e-07, 4.1502e-10, 2.3464e-12, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     3.3590e-10, 4.5952e-04, 2.4009e-07, 5.2195e-10, 2.9539e-12, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     5.2934e-10, 4.5952e-04, 2.3988e-07, 6.5619e-10, 3.7184e-12, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     8.3318e-10, 4.5952e-04, 2.3963e-07, 8.2470e-10, 4.6809e-12, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     1.3092e-09, 4.5952e-04, 2.3931e-07, 1.0361e-09, 5.8923e-12, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     2.0520e-09, 4.5952e-04, 2.3891e-07, 1.3008e-09, 7.4171e-12, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     3.2069e-09, 4.5952e-04, 2.3840e-07, 1.6320e-09, 9.3362e-12, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
     4.9947e-09, 4.5952e-04, 2.3776e-07, 2.0459e-09, 1.1751e-11, 5.7743e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7664e-04,
     7.7432e-09, 4.5952e-04, 2.3696e-07, 2.5618e-09, 1.4790e-11, 5.7742e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7664e-04,
     1.1937e-08, 4.5951e-04, 2.3596e-07, 3.2033e-09, 1.8614e-11, 5.7742e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7665e-04,
     1.8274e-08, 4.5951e-04, 2.3471e-07, 3.9989e-09, 2.3425e-11, 5.7742e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7666e-04,
     2.7726e-08, 4.5950e-04, 2.3315e-07, 4.9812e-09, 2.9476e-11, 5.7741e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7668e-04,
     4.1615e-08, 4.5949e-04, 2.3123e-07, 6.1883e-09, 3.7083e-11, 5.7740e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7670e-04,
     6.1624e-08, 4.5947e-04, 2.2885e-07, 7.6627e-09, 4.6647e-11, 5.7740e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7673e-04,
     8.9785e-08, 4.5945e-04, 2.2594e-07, 9.4510e-09, 5.8665e-11, 5.7739e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7677e-04,
     1.2820e-07, 4.5941e-04, 2.2239e-07, 1.1599e-08, 7.3758e-11, 5.7738e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7683e-04,
     1.7881e-07, 4.5937e-04, 2.1813e-07, 1.4155e-08, 9.2704e-11, 5.7736e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7690e-04,
     2.4238e-07, 4.5931e-04, 2.1303e-07, 1.7153e-08, 1.1647e-10, 5.7735e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7698e-04,
     3.1787e-07, 4.5924e-04, 2.0704e-07, 2.0620e-08, 1.4626e-10, 5.7733e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7707e-04,
     4.0141e-07, 4.5916e-04, 2.0009e-07, 2.4558e-08, 1.8358e-10, 5.7731e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7716e-04,
     4.8587e-07, 4.5908e-04, 1.9217e-07, 2.8945e-08, 2.3030e-10, 5.7729e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7725e-04,
     5.6160e-07, 4.5902e-04, 1.8334e-07, 3.3730e-08, 2.8874e-10, 5.7727e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7731e-04,
     6.1808e-07, 4.5897e-04, 1.7368e-07, 3.8832e-08, 3.6178e-10, 5.7724e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7735e-04,
     6.4739e-07, 4.5895e-04, 1.6340e-07, 4.4160e-08, 4.5301e-10, 5.7721e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7736e-04,
     6.4616e-07, 4.5896e-04, 1.5271e-07, 4.9619e-08, 5.6686e-10, 5.7719e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7734e-04,
     6.1709e-07, 4.5900e-04, 1.4189e-07, 5.5139e-08, 7.0885e-10, 5.7716e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7731e-04,
     5.6778e-07, 4.5906e-04, 1.3124e-07, 6.0702e-08, 8.8581e-10, 5.7713e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7726e-04,
     5.0793e-07, 4.5913e-04, 1.2102e-07, 6.6352e-08, 1.1062e-09, 5.7710e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7720e-04,
     4.4684e-07, 4.5920e-04, 1.1145e-07, 7.2219e-08, 1.3808e-09, 5.7707e-05],
    [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7716e-04,
     3.9146e-07, 4.5927e-04, 1.0271e-07, 7.8515e-08, 1.7228e-09, 5.7703e-05],
    [8.5370e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7712e-04,
     3.4619e-07, 4.5932e-04, 9.4879e-08, 8.5554e-08, 2.1493e-09, 5.7699e-05],
    [8.5370e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7709e-04,
     3.1323e-07, 4.5936e-04, 8.8017e-08, 9.3746e-08, 2.6819e-09, 5.7695e-05],
    [8.5369e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7708e-04,
     2.9367e-07, 4.5938e-04, 8.2121e-08, 1.0363e-07, 3.3476e-09, 5.7689e-05],
    [8.5368e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7708e-04,
     2.8858e-07, 4.5939e-04, 7.7157e-08, 1.1589e-07, 4.1815e-09, 5.7683e-05],
    [8.5367e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7709e-04,
     3.0016e-07, 4.5937e-04, 7.3078e-08, 1.3143e-07, 5.2277e-09, 5.7674e-05],
    [8.5366e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7713e-04,
     3.3284e-07, 4.5934e-04, 6.9819e-08, 1.5144e-07, 6.5418e-09, 5.7663e-05],
    [8.5365e-01, 1.4538e-01, 2.9684e-06, 1.8303e-07, 3.7719e-04,
     3.9513e-07, 4.5928e-04, 6.7313e-08, 1.7748e-07, 8.1944e-09, 5.7649e-05],
    [8.5364e-01, 1.4538e-01, 2.9684e-06, 1.8303e-07, 3.7730e-04,
     5.0250e-07, 4.5917e-04, 6.5479e-08, 2.1158e-07, 1.0273e-08, 5.7631e-05],
    [8.5364e-01, 1.4538e-01, 2.9684e-06, 1.8303e-07, 3.7748e-04,
     6.8168e-07, 4.5899e-04, 6.4227e-08, 2.5630e-07, 1.2888e-08, 5.7607e-05],
    [8.5364e-01, 1.4538e-01, 2.9684e-06, 1.8303e-07, 3.7778e-04,
     9.7748e-07, 4.5869e-04, 6.3449e-08, 3.1480e-07, 1.6173e-08, 5.7576e-05],
    [8.5365e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7828e-04,
     1.4627e-06, 4.5820e-04, 6.3029e-08, 3.9081e-07, 2.0289e-08, 5.7536e-05],
    [8.5365e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7907e-04,
     2.2498e-06, 4.5741e-04, 6.2846e-08, 4.8862e-07, 2.5424e-08, 5.7485e-05],
    [8.5365e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.8034e-04,
     3.5076e-06, 4.5615e-04, 6.2801e-08, 6.1322e-07, 3.1789e-08, 5.7420e-05],
    [8.5366e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.8233e-04,
     5.4886e-06, 4.5416e-04, 6.2830e-08, 7.7084e-07, 3.9608e-08, 5.7338e-05],
    [8.5366e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.8541e-04,
     8.5614e-06, 4.5109e-04, 6.2899e-08, 9.6929e-07, 4.9084e-08, 5.7234e-05],
    [8.5366e-01, 1.4539e-01, 2.9685e-06, 1.8304e-07, 3.9012e-04,
     1.3256e-05, 4.4639e-04, 6.2999e-08, 1.2187e-06, 6.0342e-08, 5.7104e-05],
    [8.5365e-01, 1.4539e-01, 2.9685e-06, 1.8304e-07, 3.9718e-04,
     2.0294e-05, 4.3934e-04, 6.3121e-08, 1.5317e-06, 7.3333e-08, 5.6942e-05],
    [8.5364e-01, 1.4540e-01, 2.9686e-06, 1.8304e-07, 4.0749e-04,
     3.0583e-05, 4.2905e-04, 6.3235e-08, 1.9244e-06, 8.7717e-08, 5.6740e-05],
    [8.5363e-01, 1.4540e-01, 2.9687e-06, 1.8305e-07, 4.2209e-04,
     4.5150e-05, 4.1448e-04, 6.3268e-08, 2.4165e-06, 1.0276e-07, 5.6488e-05],
    [8.5360e-01, 1.4541e-01, 2.9688e-06, 1.8306e-07, 4.4194e-04,
     6.4970e-05, 3.9467e-04, 6.3066e-08, 3.0324e-06, 1.1731e-07, 5.6176e-05],
    [8.5357e-01, 1.4541e-01, 2.9690e-06, 1.8307e-07, 4.6776e-04,
     9.0756e-05, 3.6889e-04, 6.2379e-08, 3.8024e-06, 1.2996e-07, 5.5787e-05],
    [8.5353e-01, 1.4542e-01, 2.9692e-06, 1.8308e-07, 4.9971e-04,
     1.2268e-04, 3.3700e-04, 6.0862e-08, 4.7631e-06, 1.3929e-07, 5.5306e-05],
    [8.5348e-01, 1.4544e-01, 2.9694e-06, 1.8309e-07, 5.3728e-04,
     1.6020e-04, 2.9950e-04, 5.8138e-08, 5.9590e-06, 1.4414e-07, 5.4710e-05],
    [8.5343e-01, 1.4545e-01, 2.9697e-06, 1.8311e-07, 5.7906e-04,
     2.0195e-04, 2.5780e-04, 5.3912e-08, 7.4434e-06, 1.4392e-07, 5.3973e-05],
    [8.5338e-01, 1.4546e-01, 2.9699e-06, 1.8313e-07, 6.2293e-04,
     2.4578e-04, 2.1402e-04, 4.8122e-08, 9.2794e-06, 1.3862e-07, 5.3063e-05],
    [8.5332e-01, 1.4547e-01, 2.9702e-06, 1.8314e-07, 6.6617e-04,
     2.8899e-04, 1.7087e-04, 4.1058e-08, 1.1540e-05, 1.2888e-07, 5.1943e-05],
    [8.5327e-01, 1.4549e-01, 2.9705e-06, 1.8316e-07, 7.0606e-04,
     3.2886e-04, 1.3106e-04, 3.3348e-08, 1.4306e-05, 1.1583e-07, 5.0571e-05],
    [8.5322e-01, 1.4550e-01, 2.9707e-06, 1.8317e-07, 7.4047e-04,
     3.6325e-04, 9.6725e-05, 2.5782e-08, 1.7667e-05, 1.0089e-07, 4.8902e-05],
    [8.5318e-01, 1.4551e-01, 2.9708e-06, 1.8318e-07, 7.6826e-04,
     3.9103e-04, 6.8996e-05, 1.9053e-08, 2.1713e-05, 8.5486e-08, 4.6890e-05]])

expected_radius = np.array([
    7.33195650e+09, 7.32831268e+09, 7.32467248e+09, 7.32103589e+09,
    7.31740291e+09, 7.31377352e+09, 7.31014771e+09, 7.30652549e+09,
    7.30290684e+09, 7.29929175e+09, 7.29568022e+09, 7.29207222e+09,
    7.28846775e+09, 7.28486679e+09, 7.28126932e+09, 7.27767531e+09,
    7.27408475e+09, 7.27049760e+09, 7.26691382e+09, 7.26333335e+09,
    7.25975615e+09, 7.25618213e+09, 7.25261120e+09, 7.24904325e+09,
    7.24547813e+09, 7.24191566e+09, 7.23835561e+09, 7.23479770e+09,
    7.23124157e+09, 7.22768677e+09, 7.22413274e+09, 7.22057877e+09,
    7.21702399e+09, 7.21346729e+09, 7.20990732e+09, 7.20634236e+09,
    7.20277032e+09, 7.19918862e+09, 7.19559409e+09, 7.19198290e+09,
    7.18835043e+09, 7.18469120e+09, 7.18099877e+09, 7.17726570e+09,
    7.17348355e+09, 7.16964296e+09, 7.16573383e+09, 7.16174551e+09,
    7.15766725e+09, 7.15348861e+09, 7.14920000e+09, 7.14479333e+09,
    7.14026261e+09, 7.13560460e+09, 7.13081938e+09, 7.12591084e+09,
    7.12088693e+09, 7.11575961e+09, 7.11054445e+09, 7.10525955e+09,
    7.09992411e+09, 7.09455648e+09, 7.08917228e+09, 7.08378310e+09,
    7.07839620e+09, 7.07301531e+09, 7.06764188e+09, 7.06227633e+09,
    7.05691868e+09, 7.05156887e+09, 7.04622681e+09, 7.04089240e+09,
    7.03556550e+09, 7.03024596e+09, 7.02493358e+09, 7.01962812e+09,
    7.01432929e+09, 7.00903671e+09, 7.00374992e+09, 6.99846837e+09,
    6.99319134e+09])


@pytest.mark.parametrize('call', ['-v', '--version'])
def test_command_line_version(capfd, call):
    subprocess.call(['pbay', call])
    captured = capfd.readouterr()
    assert f'Pyrat Bay version {pb.__version__}' in captured.out


def test_command_line_root(capfd):
    subprocess.call('pbay --root'.split())
    captured = capfd.readouterr()
    assert pb.constants.ROOT in captured.out


# Warm up, check when units are well or wrongly set:
def test_units_variable_not_needed(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg')
    pressure, temperature, abundances, species, radius = pb.run(cfg)


def test_units_separate(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0',
               'mpunits':'mjup'})
    pressure, temperature, abundances, species, radius = pb.run(cfg)


def test_units_in_value(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0 rjup'})
    pressure, temperature, abundances, species, radius = pb.run(cfg)


def test_units_missing(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid units 'None' for parameter mplanet." in captured.out


def test_units_invalid(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0',
               'mpunits':'nope'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid planet mass units (mpunits): nope" in captured.out


def test_units_in_value_invalid(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0 nope'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid units for value '1.0 nope' for parameter mplanet." \
        in captured.out


@pytest.mark.sort(order=1)
def test_tli_hitran_wfc3():
    pb.run(ROOT+'tests/configs/tli_hitran_1.1-1.7um_test.cfg')
    # TBD: asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_repack():
    pb.run('configs/tli_repack_test.cfg')
    # TBD: asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_tio_schwenke():
    pb.run('configs/tli_tio_schwenke_test.cfg')
    # TBD: asserts on output file


def test_pt_isothermal(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg')

    pressure, temperature, abundances, species, radius = pb.run(cfg)
    np.testing.assert_allclose(pressure, expected_pressure, rtol=1e-7)
    np.testing.assert_equal(temperature, np.tile(1500.0, 81))


def test_pt_TCEA(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_tcea.cfg')
    pressure, temperature, abundances, species, radius = pb.run(cfg)
    np.testing.assert_allclose(pressure, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temperature, expected_temperature, rtol=1e-7)


def test_atmosphere_uniform(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_uniform_test.cfg',
        reset={'atmfile':atmfile})

    press, temp, abund, species, radius = pb.run(cfg)
    np.testing.assert_allclose(press, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temp, expected_temperature, rtol=1e-7)
    q = np.tile([0.85, 0.149, 3.0e-6, 4.0e-4, 1.0e-4, 5.0e-4, 1.0e-7], (81,1))
    np.testing.assert_equal(abund, q)
    # Compare against the atmospheric file now:
    atm = io.read_atm(atmfile)
    assert atm[0] == ('bar', 'kelvin', 'volume', None)
    np.testing.assert_equal(atm[1], np.array('H2 He Na H2O CH4 CO CO2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atm[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atm[3], expected_temperature, rtol=1e-6)
    np.testing.assert_equal(atm[4], q)


def test_atmosphere_tea(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_tea_test.cfg',
        reset={'atmfile':atmfile})

    press, temp, abund, species, radius = pb.run(cfg)
    np.testing.assert_allclose(press, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temp,  expected_temperature, atol=1e-7)
    np.testing.assert_allclose(abund, expected_abundance, rtol=1e-7)
    # Compare against the atmospheric file now:
    atmf = io.read_atm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'volume', None)
    np.testing.assert_equal(atmf[1],
        np.array('H2 He Na K H2O CH4 CO CO2 NH3 HCN N2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atmf[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atmf[3], expected_temperature, rtol=5e-6)
    np.testing.assert_allclose(atmf[4], expected_abundance, rtol=1e-7)


def test_atmosphere_hydro(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_hydro_test.cfg',
        reset={'atmfile':atmfile})

    press, temp, abund, species, radius = pb.run(cfg)
    np.testing.assert_allclose(press, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temp, expected_temperature, rtol=1e-7)
    q = np.tile([0.85, 0.149, 3.0e-6, 4.0e-4, 1.0e-4, 5.0e-4, 1.0e-7], (81,1))
    np.testing.assert_equal(abund, q)
    np.testing.assert_allclose(radius, expected_radius, rtol=1e-7)
    # Compare against the atmospheric file now:
    atmf = io.read_atm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'volume', 'rjup')
    np.testing.assert_equal(atmf[1],
        np.array('H2 He Na H2O CH4 CO CO2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atmf[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atmf[3], expected_temperature, rtol=5e-6)
    np.testing.assert_equal(atmf[4], q)
    np.testing.assert_allclose(atmf[5]*pc.rjup, expected_radius, rtol=5e-5)


def test_atmosphere_hydro_default_runits(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_hydro_test.cfg',
        reset={'atmfile':atmfile, 'gplanet':'2478.7504116251885'},
        remove=['rplanet'])

    press, temp, abund, species, radius = pb.run(cfg)
    np.testing.assert_allclose(radius, expected_radius, rtol=1e-7)
    atmf = io.read_atm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'volume', 'rjup')
    np.testing.assert_allclose(atmf[5]*pc.rjup, expected_radius, rtol=5e-5)



# See tests/test_spectrum.py for spectrum tests


def test_spectrum_emission(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'rt_path':'emission', 'cpars':'-0.5'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    # TBD: implement asserts


@pytest.mark.sort(order=10)
def test_opacity_pbay(capfd):
    pyrat = pb.run(ROOT+'tests/configs/opacity_test.cfg')
    captured = capfd.readouterr()
    assert "Extinction-coefficient table written to file:" in captured.out
    assert "exttable_test_300-3000K_1.1-1.7um.npz'." in captured.out
    assert "exttable_test_300-3000K_1.1-1.7um.npz" in os.listdir('outputs/')


@pytest.mark.skip
def test_mcmc():
    pyrat = pb.run(ROOT+'tests/configs/mcmc_transmission_test.cfg')
