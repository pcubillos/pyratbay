import os
import sys
import subprocess
import pytest

import numpy as np

from conftest import make_config

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay as pb
import pyratbay.constants  as pc
import pyratbay.atmosphere as pa

os.chdir(ROOT+'tests')

expected_pressure    = np.logspace(0, 8, 81)
expected_temperature = np.array(
      [1047.04535531, 1047.04636729, 1047.04764749, 1047.04926694,
       1047.05131548, 1047.05390677, 1047.05718449, 1047.06133039,
       1047.06657429, 1047.07320679, 1047.08159536, 1047.09220465,
       1047.10562211, 1047.12259045, 1047.1440486 , 1047.17118342,
       1047.20549505, 1047.24887932, 1047.30373179, 1047.37307893,
       1047.46074334, 1047.57155184, 1047.7115971 , 1047.88856622,
       1048.11215262, 1048.39457122, 1048.75120098, 1049.2013834 ,
       1049.76941034, 1050.48573869, 1051.38847286, 1052.52515618,
       1053.95490802, 1055.75092978, 1058.00337628, 1060.82254126,
       1064.34222961, 1068.72307524, 1074.15540631, 1080.8610604 ,
       1089.09332826, 1099.13399769, 1111.28635212, 1125.86305112,
       1143.16818192, 1163.4734689 , 1186.98959518, 1213.83461222,
       1244.00218101, 1277.3326441 , 1313.48964651, 1351.94449886,
       1391.97022283, 1432.64772742, 1472.88802251, 1511.47646601,
       1547.14676037, 1578.69185734, 1605.11307825, 1625.79393924,
       1640.65971886, 1650.25479478, 1655.66159741, 1658.23445306,
       1659.23535547, 1659.55574759, 1659.66292568, 1659.73229419,
       1659.81020803, 1659.90748667, 1660.02989399, 1660.1839565 ,
       1660.37784873, 1660.62184804, 1660.92887213, 1661.31515062,
       1661.80106362, 1662.4121864 , 1663.18058734, 1664.14643497,
       1665.35997887])

expected_abundance = np.array(
    [[8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      8.5798e-15, 4.5952e-04, 2.4071e-07, 2.6266e-12, 1.4806e-14, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      1.3597e-14, 4.5952e-04, 2.4071e-07, 3.3066e-12, 1.8640e-14, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      2.1552e-14, 4.5952e-04, 2.4071e-07, 4.1629e-12, 2.3467e-14, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      3.4158e-14, 4.5952e-04, 2.4071e-07, 5.2409e-12, 2.9543e-14, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      5.4136e-14, 4.5952e-04, 2.4071e-07, 6.5978e-12, 3.7192e-14, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      8.5800e-14, 4.5952e-04, 2.4071e-07, 8.3062e-12, 4.6823e-14, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      1.3595e-13, 4.5952e-04, 2.4070e-07, 1.0456e-11, 5.8946e-14, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      2.1547e-13, 4.5952e-04, 2.4070e-07, 1.3164e-11, 7.4209e-14, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      3.4141e-13, 4.5952e-04, 2.4069e-07, 1.6571e-11, 9.3423e-14, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      5.4109e-13, 4.5952e-04, 2.4069e-07, 2.0862e-11, 1.1761e-13, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      8.5735e-13, 4.5952e-04, 2.4068e-07, 2.6262e-11, 1.4806e-13, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      1.3584e-12, 4.5952e-04, 2.4068e-07, 3.3059e-11, 1.8640e-13, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      2.1520e-12, 4.5952e-04, 2.4066e-07, 4.1614e-11, 2.3467e-13, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      3.4100e-12, 4.5952e-04, 2.4065e-07, 5.2387e-11, 2.9543e-13, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      5.4016e-12, 4.5952e-04, 2.4063e-07, 6.5943e-11, 3.7192e-13, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      8.5546e-12, 4.5952e-04, 2.4060e-07, 8.3002e-11, 4.6821e-13, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      1.3545e-11, 4.5952e-04, 2.4057e-07, 1.0447e-10, 5.8944e-13, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      2.1446e-11, 4.5952e-04, 2.4053e-07, 1.3149e-10, 7.4205e-13, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      3.3947e-11, 4.5952e-04, 2.4048e-07, 1.6548e-10, 9.3418e-13, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      5.3710e-11, 4.5952e-04, 2.4042e-07, 2.0824e-10, 1.1760e-12, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      8.4935e-11, 4.5952e-04, 2.4034e-07, 2.6202e-10, 1.4805e-12, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      1.3424e-10, 4.5952e-04, 2.4024e-07, 3.2963e-10, 1.8638e-12, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      2.1204e-10, 4.5952e-04, 2.4011e-07, 4.1465e-10, 2.3463e-12, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      3.3458e-10, 4.5952e-04, 2.3995e-07, 5.2145e-10, 2.9537e-12, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      5.2739e-10, 4.5952e-04, 2.3975e-07, 6.5559e-10, 3.7182e-12, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      8.3011e-10, 4.5952e-04, 2.3950e-07, 8.2396e-10, 4.6807e-12, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      1.3040e-09, 4.5952e-04, 2.3917e-07, 1.0351e-09, 5.8921e-12, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      2.0440e-09, 4.5952e-04, 2.3877e-07, 1.2995e-09, 7.4168e-12, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      3.1944e-09, 4.5952e-04, 2.3825e-07, 1.6305e-09, 9.3357e-12, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7663e-04,
      4.9739e-09, 4.5952e-04, 2.3761e-07, 2.0438e-09, 1.1751e-11, 5.7743e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7664e-04,
      7.7110e-09, 4.5952e-04, 2.3681e-07, 2.5592e-09, 1.4790e-11, 5.7742e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7664e-04,
      1.1884e-08, 4.5951e-04, 2.3580e-07, 3.1999e-09, 1.8613e-11, 5.7742e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7665e-04,
      1.8194e-08, 4.5951e-04, 2.3455e-07, 3.9946e-09, 2.3424e-11, 5.7742e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7666e-04,
      2.7599e-08, 4.5950e-04, 2.3299e-07, 4.9755e-09, 2.9474e-11, 5.7741e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7668e-04,
      4.1414e-08, 4.5949e-04, 2.3106e-07, 6.1809e-09, 3.7081e-11, 5.7740e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7670e-04,
      6.1313e-08, 4.5947e-04, 2.2867e-07, 7.6532e-09, 4.6644e-11, 5.7740e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7673e-04,
      8.9292e-08, 4.5945e-04, 2.2575e-07, 9.4383e-09, 5.8661e-11, 5.7739e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7677e-04,
      1.2747e-07, 4.5941e-04, 2.2220e-07, 1.1583e-08, 7.3753e-11, 5.7738e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7683e-04,
      1.7764e-07, 4.5937e-04, 2.1791e-07, 1.4132e-08, 9.2697e-11, 5.7736e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7690e-04,
      2.4070e-07, 4.5931e-04, 2.1281e-07, 1.7124e-08, 1.1646e-10, 5.7735e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7698e-04,
      3.1549e-07, 4.5924e-04, 2.0680e-07, 2.0582e-08, 1.4625e-10, 5.7733e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7707e-04,
      3.9810e-07, 4.5916e-04, 1.9984e-07, 2.4508e-08, 1.8357e-10, 5.7731e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7716e-04,
      4.8141e-07, 4.5909e-04, 1.9191e-07, 2.8879e-08, 2.3028e-10, 5.7729e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7724e-04,
      5.5599e-07, 4.5902e-04, 1.8306e-07, 3.3647e-08, 2.8871e-10, 5.7727e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7731e-04,
      6.1132e-07, 4.5898e-04, 1.7340e-07, 3.8727e-08, 3.6174e-10, 5.7724e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7735e-04,
      6.3978e-07, 4.5896e-04, 1.6311e-07, 4.4032e-08, 4.5295e-10, 5.7721e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7735e-04,
      6.3799e-07, 4.5897e-04, 1.5242e-07, 4.9463e-08, 5.6679e-10, 5.7719e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7734e-04,
      6.0895e-07, 4.5901e-04, 1.4162e-07, 5.4958e-08, 7.0875e-10, 5.7716e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7730e-04,
      5.5994e-07, 4.5907e-04, 1.3098e-07, 6.0494e-08, 8.8566e-10, 5.7713e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7725e-04,
      5.0076e-07, 4.5914e-04, 1.2078e-07, 6.6119e-08, 1.1060e-09, 5.7710e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7720e-04,
      4.4045e-07, 4.5921e-04, 1.1123e-07, 7.1962e-08, 1.3805e-09, 5.7707e-05],
     [8.5371e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7715e-04,
      3.8593e-07, 4.5927e-04, 1.0251e-07, 7.8239e-08, 1.7224e-09, 5.7703e-05],
     [8.5370e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7711e-04,
      3.4142e-07, 4.5932e-04, 9.4701e-08, 8.5261e-08, 2.1489e-09, 5.7700e-05],
     [8.5370e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.7709e-04,
      3.0910e-07, 4.5936e-04, 8.7862e-08, 9.3439e-08, 2.6813e-09, 5.7695e-05],
     [8.5369e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7707e-04,
      2.9004e-07, 4.5938e-04, 8.1987e-08, 1.0331e-07, 3.3470e-09, 5.7690e-05],
     [8.5368e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7707e-04,
      2.8533e-07, 4.5939e-04, 7.7044e-08, 1.1556e-07, 4.1808e-09, 5.7683e-05],
     [8.5367e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7709e-04,
      2.9714e-07, 4.5938e-04, 7.2982e-08, 1.3110e-07, 5.2269e-09, 5.7674e-05],
     [8.5366e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7712e-04,
      3.2995e-07, 4.5934e-04, 6.9742e-08, 1.5112e-07, 6.5410e-09, 5.7663e-05],
     [8.5365e-01, 1.4538e-01, 2.9684e-06, 1.8303e-07, 3.7719e-04,
      3.9227e-07, 4.5928e-04, 6.7251e-08, 1.7716e-07, 8.1935e-09, 5.7649e-05],
     [8.5364e-01, 1.4538e-01, 2.9684e-06, 1.8303e-07, 3.7730e-04,
      4.9953e-07, 4.5917e-04, 6.5430e-08, 2.1127e-07, 1.0273e-08, 5.7631e-05],
     [8.5364e-01, 1.4538e-01, 2.9684e-06, 1.8303e-07, 3.7748e-04,
      6.7842e-07, 4.5899e-04, 6.4188e-08, 2.5600e-07, 1.2888e-08, 5.7607e-05],
     [8.5364e-01, 1.4538e-01, 2.9684e-06, 1.8303e-07, 3.7778e-04,
      9.7385e-07, 4.5869e-04, 6.3420e-08, 3.1451e-07, 1.6172e-08, 5.7576e-05],
     [8.5365e-01, 1.4538e-01, 2.9684e-06, 1.8303e-07, 3.7827e-04,
      1.4582e-06, 4.5821e-04, 6.3004e-08, 3.9051e-07, 2.0288e-08, 5.7536e-05],
     [8.5365e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.7906e-04,
      2.2438e-06, 4.5742e-04, 6.2825e-08, 4.8829e-07, 2.5423e-08, 5.7485e-05],
     [8.5365e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.8033e-04,
      3.4990e-06, 4.5616e-04, 6.2781e-08, 6.1284e-07, 3.1789e-08, 5.7420e-05],
     [8.5366e-01, 1.4539e-01, 2.9684e-06, 1.8303e-07, 3.8231e-04,
      5.4753e-06, 4.5418e-04, 6.2810e-08, 7.7036e-07, 3.9609e-08, 5.7338e-05],
     [8.5366e-01, 1.4539e-01, 2.9685e-06, 1.8303e-07, 3.8539e-04,
      8.5416e-06, 4.5111e-04, 6.2880e-08, 9.6871e-07, 4.9087e-08, 5.7234e-05],
     [8.5366e-01, 1.4539e-01, 2.9685e-06, 1.8304e-07, 3.9009e-04,
      1.3226e-05, 4.4642e-04, 6.2980e-08, 1.2179e-06, 6.0348e-08, 5.7105e-05],
     [8.5365e-01, 1.4539e-01, 2.9685e-06, 1.8304e-07, 3.9714e-04,
      2.0249e-05, 4.3939e-04, 6.3101e-08, 1.5308e-06, 7.3346e-08, 5.6943e-05],
     [8.5364e-01, 1.4540e-01, 2.9686e-06, 1.8304e-07, 4.0743e-04,
      3.0519e-05, 4.2911e-04, 6.3216e-08, 1.9232e-06, 8.7740e-08, 5.6741e-05],
     [8.5363e-01, 1.4540e-01, 2.9687e-06, 1.8305e-07, 4.2200e-04,
      4.5060e-05, 4.1457e-04, 6.3249e-08, 2.4150e-06, 1.0280e-07, 5.6489e-05],
     [8.5360e-01, 1.4541e-01, 2.9688e-06, 1.8306e-07, 4.4182e-04,
      6.4851e-05, 3.9479e-04, 6.3049e-08, 3.0306e-06, 1.1737e-07, 5.6176e-05],
     [8.5357e-01, 1.4541e-01, 2.9690e-06, 1.8307e-07, 4.6760e-04,
      9.0598e-05, 3.6905e-04, 6.2365e-08, 3.8001e-06, 1.3006e-07, 5.5788e-05],
     [8.5353e-01, 1.4542e-01, 2.9692e-06, 1.8308e-07, 4.9953e-04,
      1.2249e-04, 3.3718e-04, 6.0855e-08, 4.7603e-06, 1.3941e-07, 5.5308e-05],
     [8.5348e-01, 1.4544e-01, 2.9694e-06, 1.8309e-07, 5.3706e-04,
      1.5999e-04, 2.9972e-04, 5.8139e-08, 5.9554e-06, 1.4430e-07, 5.4712e-05],
     [8.5343e-01, 1.4545e-01, 2.9697e-06, 1.8311e-07, 5.7883e-04,
      2.0172e-04, 2.5803e-04, 5.3922e-08, 7.4389e-06, 1.4411e-07, 5.3975e-05],
     [8.5338e-01, 1.4546e-01, 2.9699e-06, 1.8312e-07, 6.2269e-04,
      2.4554e-04, 2.1426e-04, 4.8142e-08, 9.2739e-06, 1.3883e-07, 5.3066e-05],
     [8.5332e-01, 1.4547e-01, 2.9702e-06, 1.8314e-07, 6.6594e-04,
      2.8876e-04, 1.7110e-04, 4.1086e-08, 1.1533e-05, 1.2910e-07, 5.1946e-05],
     [8.5327e-01, 1.4549e-01, 2.9705e-06, 1.8316e-07, 7.0586e-04,
      3.2866e-04, 1.3126e-04, 3.3380e-08, 1.4298e-05, 1.1604e-07, 5.0575e-05],
     [8.5322e-01, 1.4550e-01, 2.9707e-06, 1.8317e-07, 7.4029e-04,
      3.6307e-04, 9.6904e-05, 2.5815e-08, 1.7657e-05, 1.0110e-07, 4.8907e-05],
     [8.5318e-01, 1.4551e-01, 2.9708e-06, 1.8318e-07, 7.6812e-04,
      3.9089e-04, 6.9137e-05, 1.9083e-08, 2.1701e-05, 8.5679e-08, 4.6897e-05]])


# Warm up, check when units are well or wrongly set:
def test_units_variable_not_needed(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg')
    pressure, temperature = pb.run(cfg)


def test_units_separate(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={'mplanet':'1.0',
               'mpunits':'mjup'})
    pressure, temperature = pb.run(cfg)


def test_units_in_value(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={'mplanet':'1.0 rjup'})
    pressure, temperature = pb.run(cfg)


def test_units_missing(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={'mplanet':'1.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid units 'None' for parameter mplanet." in captured.out


def test_units_invalid(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={'mplanet':'1.0',
               'mpunits':'nope'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid planet mass units (mpunits): nope" in captured.out


def test_units_in_value_invalid(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={'mplanet':'1.0 nope'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid units for value '1.0 nope' of parameter mplanet." \
        in captured.out


@pytest.mark.sort(order=1)
def test_tli_hitran_wfc3():
    pb.run(ROOT+'tests/tli_hitran_1.1-1.7um_test.cfg')
    # TBD: asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_pands():
    pb.run('tli_pands_test.cfg')
    # TBD: asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_repack():
    pb.run('tli_repack_test.cfg')
    # TBD: asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_tio_schwenke():
    pb.run('tli_tio_schwenke_test.cfg')
    # TBD: asserts on output file


def test_pt_isothermal(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg')

    pressure, temperature = pb.run(cfg)
    np.testing.assert_allclose(pressure, expected_pressure, rtol=1e-7)
    np.testing.assert_equal(temperature, np.tile(1500.0, 81))


def test_pt_TCEA(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/pt_tcea.cfg')

    pressure, temperature = pb.run(cfg)
    np.testing.assert_allclose(pressure, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temperature, expected_temperature, atol=1e-10)


def test_atmosphere_uniform(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path, ROOT+'tests/atmosphere_uniform_test.cfg',
        reset={'atmfile':atmfile})

    atm = pb.run(cfg)
    np.testing.assert_allclose(atm[0], expected_pressure,    rtol=1e-7)
    np.testing.assert_allclose(atm[1], expected_temperature, atol=1e-10)
    q = np.tile([0.85, 0.149, 3.0e-6, 4.0e-4, 1.0e-4, 5.0e-4, 1.0e-7], (81,1))
    np.testing.assert_equal(atm[2], q)
    # Compare against the atmospheric file now:
    atm = pa.readatm(atmfile)
    assert atm[0] == ('bar', 'kelvin', 'number', None)
    np.testing.assert_equal(atm[1], np.array('H2 He Na H2O CH4 CO CO2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atm[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atm[3], expected_temperature,     atol=5e-4)
    np.testing.assert_equal(atm[4], q)


def test_atmosphere_tea(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path, ROOT+'tests/atmosphere_tea_test.cfg',
        reset={'atmfile':atmfile})

    atm = pb.run(cfg)
    np.testing.assert_allclose(atm[0], expected_pressure,    rtol=1e-7)
    np.testing.assert_allclose(atm[1], expected_temperature, atol=1e-10)
    np.testing.assert_allclose(atm[2], expected_abundance,   rtol=1e-7)
    # Compare against the atmospheric file now:
    atmf = pa.readatm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'number', None)
    np.testing.assert_equal(atmf[1],
        np.array('H2 He Na K H2O CH4 CO CO2 NH3 HCN N2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atmf[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atmf[3], expected_temperature,     atol=5e-3)
    np.testing.assert_allclose(atmf[4], expected_abundance,       rtol=1e-7)


# See tests/test_spectrum.py for spectrum tests


def test_spectrum_emission(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'path':'eclipse', 'cpars':'-0.5'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    # TBD: implement asserts


@pytest.mark.sort(order=10)
def test_opacity(capfd):
    pyrat = pb.run(ROOT+'tests/opacity_test.cfg')
    captured = capfd.readouterr()
    assert "Extinction-coefficient table written to file:" in captured.out
    assert "exttable_test_300-3000K_1.1-1.7um.dat'." in captured.out
    assert "exttable_test_300-3000K_1.1-1.7um.dat" in os.listdir('.')


@pytest.mark.skip
def test_mcmc():
    pyrat = pb.run(ROOT+'tests/mcmc_transmission_test.cfg')
