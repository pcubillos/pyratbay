# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest
import re
import shutil

import numpy as np
import pyratbay as pb
import pyratbay.atmosphere as pa
import pyratbay.constants as pc
from pyratbay.constants import ROOT

from conftest import make_config

os.chdir(ROOT + 'tests')


pmin, pmax = 0, 8
# Look at tests/configs/atmosphere_jupiter_read.cfg
read_t = 1000.0
read_nlayers = 51
# Look at tests/configs/atmosphere_jupiter_calc.cfg
calc_t = 1500.0
calc_nlayers = 71


#expected_vmr = np.load(
#    f'{ROOT}/tests/expected/expected_tea_profile.npz')['arr_0']

expected_radius = np.array([
    7.40887177e+09, 7.40349361e+09, 7.39812325e+09, 7.39276067e+09,
    7.38740586e+09, 7.38205881e+09, 7.37671949e+09, 7.37138789e+09,
    7.36606399e+09, 7.36074777e+09, 7.35543922e+09, 7.35013832e+09,
    7.34484506e+09, 7.33955942e+09, 7.33428138e+09, 7.32901092e+09,
    7.32374804e+09, 7.31849270e+09, 7.31324491e+09, 7.30800463e+09,
    7.30277186e+09, 7.29754658e+09, 7.29232877e+09, 7.28711841e+09,
    7.28191550e+09, 7.27672001e+09, 7.27153193e+09, 7.26635124e+09,
    7.26117793e+09, 7.25601197e+09, 7.25085337e+09, 7.24570209e+09,
    7.24055813e+09, 7.23542147e+09, 7.23029209e+09, 7.22516998e+09,
    7.22005511e+09, 7.21494749e+09, 7.20984709e+09, 7.20475389e+09,
    7.19966789e+09, 7.19458906e+09, 7.18951739e+09, 7.18445286e+09,
    7.17939547e+09, 7.17434519e+09, 7.16930201e+09, 7.16426592e+09,
    7.15923690e+09, 7.15421493e+09, 7.14920000e+09, 7.14419210e+09,
    7.13919121e+09, 7.13419731e+09, 7.12921040e+09, 7.12423046e+09,
    7.11925746e+09, 7.11429141e+09, 7.10933227e+09, 7.10438005e+09,
    7.09943472e+09, 7.09449627e+09, 7.08956469e+09, 7.08463996e+09,
    7.07972207e+09, 7.07481099e+09, 7.06990673e+09, 7.06500926e+09,
    7.06011857e+09, 7.05523465e+09, 7.05035748e+09, 7.04548705e+09,
    7.04062335e+09, 7.03576635e+09, 7.03091605e+09, 7.02607244e+09,
    7.02123549e+09, 7.01640520e+09, 7.01158155e+09, 7.00676452e+09,
    7.00195411e+09,
])

expected_calc_radius = pc.rjup * np.array([
 1.5831851 , 1.58117967, 1.57917931, 1.57718401, 1.57519375, 1.5732085 ,
 1.57122825, 1.56925298, 1.56728267, 1.5653173 , 1.56335685, 1.56140131,
 1.55945065, 1.55750486, 1.55556392, 1.55362782, 1.55169652, 1.54977002,
 1.5478483 , 1.54593134, 1.54401912, 1.54211163, 1.54020884, 1.53831075,
 1.53641732, 1.53452855, 1.53264442, 1.53076491, 1.52889001, 1.52701969,
 1.52515394, 1.52329275, 1.52143609, 1.51958396, 1.51773632, 1.51589318,
 1.5140545 , 1.51222028, 1.5103905 , 1.50856515, 1.5067442 , 1.50492764,
 1.50311545, 1.50130763, 1.49950414, 1.49770499, 1.49591015, 1.4941196 ,
 1.49233333, 1.49055134, 1.48877359, 1.48700008, 1.48523078, 1.4834657 ,
 1.4817048 , 1.47994808, 1.47819552, 1.4764471 , 1.47470282, 1.47296265,
 1.47122659, 1.46949461, 1.46776671, 1.46604286, 1.46432306, 1.46260729,
 1.46089554, 1.45918778, 1.45748402, 1.45578423, 1.4540884 ,
])
# read p, t, vmr, r:
expected_read_radius = pc.rjup * np.array([
 1.023869 , 1.023087 , 1.022306 , 1.021526 , 1.020747 , 1.01997  , 1.019193 ,
 1.018418 , 1.017644 , 1.016871 , 1.0161   , 1.015329 , 1.01456  , 1.013792 ,
 1.013025 , 1.012259 , 1.011494 , 1.010731 , 1.009969 , 1.009207 , 1.008447 ,
 1.007688 , 1.006931 , 1.006174 , 1.005419 , 1.004664 , 1.003911 , 1.003159 ,
 1.002408 , 1.001658 , 1.000909 , 1.000162 , 0.9994153, 0.9986699, 0.9979257,
 0.9971825, 0.9964405, 0.9956995, 0.9949597, 0.9942209, 0.9934833, 0.9927467,
 0.9920113, 0.9912769, 0.9905436, 0.9898114, 0.9890803, 0.9883502, 0.9876213,
 0.9868934, 0.9861666,
])
# read p, calc t vmr r:
expected_read_p_radius = pc.rjup * np.array([
 1.58314175, 1.58033572, 1.57753963, 1.57475341, 1.57197701, 1.56921039,
 1.56645349, 1.56370625, 1.56096864, 1.5582406 , 1.55552207, 1.55281302,
 1.55011338, 1.54742312, 1.54474217, 1.5420705 , 1.53940806, 1.53675479,
 1.53411065, 1.5314756 , 1.52884958, 1.52623256, 1.52362447, 1.52102529,
 1.51843496, 1.51585344, 1.51328067, 1.51071663, 1.50816126, 1.50561453,
 1.50307638, 1.50054677, 1.49802566, 1.49551301, 1.49300877, 1.49051291,
 1.48802538, 1.48554614, 1.48307514, 1.48061235, 1.47815773, 1.47571123,
 1.47327282, 1.47084245, 1.46842009, 1.46600569, 1.46359922, 1.46120064,
 1.45880991, 1.45642698, 1.45405183,
])
# read t, calc p vmr r:
expected_read_t_radius = pc.rjup * np.array([
 1.55445019, 1.55316076, 1.55187348, 1.55058832, 1.54930529, 1.54802438,
 1.54674559, 1.54546891, 1.54419433, 1.54292186, 1.54165148, 1.54038319,
 1.53911698, 1.53785286, 1.53659081, 1.53533083, 1.53407291, 1.53281706,
 1.53156326, 1.53031151, 1.5290618 , 1.52781413, 1.5265685 , 1.52532489,
 1.52408332, 1.52284376, 1.52160621, 1.52037068, 1.51913715, 1.51790562,
 1.51667608, 1.51544854, 1.51422298, 1.5129994 , 1.51177779, 1.51055816,
 1.5093405 , 1.50812479, 1.50691105, 1.50569925, 1.5044894 , 1.5032815 ,
 1.50207553, 1.5008715 , 1.49966939, 1.49846921, 1.49727095, 1.4960746 ,
 1.49488017, 1.49368764, 1.49249701, 1.49130827, 1.49012143, 1.48893648,
 1.48775341, 1.48657222, 1.4853929 , 1.48421545, 1.48303987, 1.48186615,
 1.48069428, 1.47952427, 1.4783561 , 1.47718978, 1.4760253 , 1.47486265,
 1.47370183, 1.47254284, 1.47138567, 1.47023031, 1.46907677,
])
# read pt, calc vmr r:
expected_read_pt_radius = pc.rjup * np.array([
 1.55442233, 1.5526178 , 1.55081745, 1.54902127, 1.54722925, 1.54544137,
 1.54365761, 1.54187797, 1.54010243, 1.53833097, 1.53656359, 1.53480026,
 1.53304097, 1.53128571, 1.52953446, 1.52778722, 1.52604396, 1.52430468,
 1.52256936, 1.52083798, 1.51911054, 1.51738701, 1.5156674 , 1.51395167,
 1.51223983, 1.51053185, 1.50882773, 1.50712745, 1.50543099, 1.50373835,
 1.50204952, 1.50036447, 1.4986832 , 1.49700569, 1.49533193, 1.49366191,
 1.49199562, 1.49033304, 1.48867416, 1.48701897, 1.48536746, 1.48371961,
 1.48207542, 1.48043486, 1.47879794, 1.47716462, 1.47553491, 1.4739088 ,
 1.47228626, 1.4706673 , 1.46905189,
])
# read p,t,vmr, calc r:
expected_read_ptq_radius = pc.rjup * np.array([
 1.55435409, 1.55255189, 1.55075387, 1.54896001, 1.5471703 , 1.54538471,
 1.54360325, 1.54182588, 1.5400526 , 1.5382834 , 1.53651826, 1.53475717,
 1.5330001 , 1.53124706, 1.52949802, 1.52775297, 1.5260119 , 1.52427479,
 1.52254163, 1.52081241, 1.51908712, 1.51736573, 1.51564824, 1.51393463,
 1.5122249 , 1.51051902, 1.50881698, 1.50711878, 1.5054244 , 1.50373382,
 1.50204703, 1.50036403, 1.49868479, 1.4970093 , 1.49533756, 1.49366955,
 1.49200525, 1.49034466, 1.48868777, 1.48703455, 1.485385  , 1.4837391 ,
 1.48209685, 1.48045823, 1.47882323, 1.47719184, 1.47556404, 1.47393983,
 1.47231918, 1.4707021 , 1.46908857,
])
# calc p, interp t,vmr, r:
expected_interp_tq_radius = pc.rjup * np.array([
 1.023869  , 1.02331043, 1.02275229, 1.02219457, 1.02163743, 1.02108086,
 1.020525  , 1.01997   , 1.019415  , 1.01886086, 1.01830743, 1.01775457,
 1.01720229, 1.01665071, 1.0161    , 1.01554929, 1.01499943, 1.01445029,
 1.01390171, 1.01335371, 1.01280614, 1.012259  , 1.01171257, 1.011167  ,
 1.01062214, 1.01007786, 1.00953357, 1.00898986, 1.008447  , 1.00790486,
 1.00736357, 1.00682286, 1.00628214, 1.00574257, 1.00520329, 1.004664  ,
 1.00412614, 1.00358871, 1.00305171, 1.00251529, 1.00197943, 1.001444  ,
 1.000909  , 1.00037543, 0.99984199, 0.99930881, 0.99877639, 0.99824464,
 0.99771336, 0.9971825 , 0.9966525 , 0.99612293, 0.99559381, 0.99506539,
 0.99453753, 0.99401016, 0.9934833 , 0.99295716, 0.99243153, 0.99190639,
 0.99138181, 0.99085787, 0.9903344 , 0.9898114 , 0.98928919, 0.9887674 ,
 0.98824607, 0.98772543, 0.98720536, 0.98668574, 0.9861666 ,
])
# calc p, interp t,vmr, calc r:
expected_interp_tq_calc_radius = pc.rjup * np.array([
 1.55438191, 1.55309415, 1.55180853, 1.55052504, 1.54924366, 1.5479644 ,
 1.54668726, 1.54541222, 1.54413927, 1.54286843, 1.54159967, 1.540333  ,
 1.53906841, 1.5378059 , 1.53654545, 1.53528707, 1.53403075, 1.53277648,
 1.53152426, 1.53027408, 1.52902595, 1.52777985, 1.52653578, 1.52529373,
 1.5240537 , 1.52281569, 1.52157969, 1.52034569, 1.51911369, 1.51788369,
 1.51665568, 1.51542965, 1.5142056 , 1.51298353, 1.51176343, 1.5105453 ,
 1.50932912, 1.50811491, 1.50690264, 1.50569233, 1.50448395, 1.50327752,
 1.50207302, 1.50087044, 1.49966979, 1.49847106, 1.49727425, 1.49607934,
 1.49488634, 1.49369524, 1.49250604, 1.49131873, 1.49013331, 1.48894977,
 1.48776811, 1.48658832, 1.48541041, 1.48423435, 1.48306016, 1.48188783,
 1.48071735, 1.47954871, 1.47838192, 1.47721697, 1.47605385, 1.47489256,
 1.4737331 , 1.47257546, 1.47141963, 1.47026562, 1.46911342,
])


@pytest.fixture
def reset_jupiter():
    shutil.copy(
        'inputs/jupiter_isothermal_uniform_vmr_read.atm',
        'inputs/jupiter_isothermal_uniform_vmr.atm',
    )


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Any input:
@pytest.mark.parametrize('atm_input', ('none', 'pt'))
def test_run_atmosphere_calc_pt(tmp_path, atm_input, reset_jupiter):
    reset = {}
    remove = ['chemistry', 'radmodel']

    if atm_input == 'none':
        remove += ['atmfile']
    elif atm_input == 'pt':
        reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)
    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', calc_nlayers)
    expected_temperature = np.tile(calc_t, calc_nlayers)

    np.testing.assert_allclose(atm_model.press, expected_pressure)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    assert atm_model.species is None
    assert atm_model.vmr is None
    assert atm_model.radius is None


@pytest.mark.parametrize('atm_input', ('none', 'pt', 'atm'))
def test_run_atmosphere_calc_ptq(tmp_path, atm_input, reset_jupiter):
    reset = {}
    remove = ['radmodel']

    if atm_input == 'none':
        remove += ['atmfile']
    elif atm_input == 'pt':
        reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'
    elif atm_input == 'atm':
        reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)
    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', calc_nlayers)
    expected_temperature = np.tile(calc_t, calc_nlayers)
    expected_species = ['H2', 'He', 'H2O']
    expected_vmr = np.array([0.85, 0.149, 1e-4])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr, rtol=1e-7)
    assert atm_model.radius is None


@pytest.mark.parametrize('atm_input', ('none', 'pt', 'atm'))
def test_run_atmosphere_calc_ptqr(tmp_path, atm_input, reset_jupiter):
    reset = {}
    remove = []
    if atm_input == 'none':
        remove += ['atmfile']
    elif atm_input == 'pt':
        reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'
    elif atm_input == 'atm':
        reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)
    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', calc_nlayers)
    expected_temperature = np.tile(calc_t, calc_nlayers)
    expected_species = ['H2', 'He', 'H2O']
    expected_vmr = np.array([0.85, 0.149, 1e-4])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)
    np.testing.assert_allclose(atm_model.radius, expected_calc_radius)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# From pt or atm file:
@pytest.mark.parametrize('atm_input', ('pt', 'atm'))
def test_run_atmosphere_read_p_calc_tq(tmp_path, atm_input, reset_jupiter):
    reset = {}
    remove = ['nlayers', 'radmodel']
    if atm_input == 'pt':
        reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'
    elif atm_input == 'atm':
        reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', read_nlayers)
    expected_temperature = np.tile(calc_t, read_nlayers)
    expected_species = ['H2', 'He', 'H2O']
    expected_vmr = np.array([0.85, 0.149, 1e-4])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)
    assert atm_model.radius is None


@pytest.mark.parametrize('atm_input', ('pt', 'atm'))
def test_run_atmosphere_read_pt_calc_qr(tmp_path, atm_input, reset_jupiter):
    reset = {}
    remove = ['nlayers', 'tmodel']
    if atm_input == 'pt':
        reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'
    elif atm_input == 'atm':
        reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', read_nlayers)
    expected_temperature = np.tile(read_t, read_nlayers)
    expected_species = ['H2', 'He', 'H2O']
    expected_vmr = np.array([0.85, 0.149, 1e-4])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)
    np.testing.assert_allclose(atm_model.radius, expected_read_pt_radius)


@pytest.mark.parametrize('atm_input', ('pt', 'atm'))
def test_run_atmosphere_read_p_calc_tqr(tmp_path, atm_input, reset_jupiter):
    reset = {}
    remove = ['nlayers']
    if atm_input == 'pt':
        reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'
    elif atm_input == 'atm':
        reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', read_nlayers)
    expected_temperature = np.tile(calc_t, read_nlayers)
    expected_species = ['H2', 'He', 'H2O']
    expected_vmr = np.array([0.85, 0.149, 1e-4])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)
    np.testing.assert_allclose(atm_model.radius, expected_read_p_radius)


@pytest.mark.parametrize('atm_input', ('pt', 'atm'))
def test_run_atmosphere_calc_pq_interp_t(tmp_path, atm_input, reset_jupiter):
    reset = {}
    remove = ['tmodel', 'radmodel']
    if atm_input == 'pt':
        reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'
    elif atm_input == 'atm':
        reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', calc_nlayers)
    expected_temperature = np.tile(read_t, calc_nlayers)
    atm_model = pb.run(cfg)
    expected_species = ['H2', 'He', 'H2O']
    expected_vmr = np.array([0.85, 0.149, 1e-4])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)
    assert atm_model.radius is None


@pytest.mark.parametrize('atm_input', ('pt', 'atm'))
def test_run_atmosphere_calc_pqr_interp_t(tmp_path, atm_input, reset_jupiter):
    reset = {}
    remove = ['tmodel']
    if atm_input == 'pt':
        reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'
    elif atm_input == 'atm':
        reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)
    expected_species = ['H2', 'He', 'H2O']
    expected_vmr = np.array([0.85, 0.149, 1e-4])
    output_vmr = atm_model.vmr[0]

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', calc_nlayers)
    expected_temperature = np.tile(read_t, calc_nlayers)

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)
    np.testing.assert_allclose(atm_model.radius, expected_read_t_radius)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# From pt file only:
def test_run_atmosphere_read_pt_from_pt(tmp_path, reset_jupiter):
    reset = {}
    remove = ['nlayers', 'tmodel', 'chemistry', 'radmodel']
    reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', read_nlayers)
    expected_temperature = np.tile(read_t, read_nlayers)

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    assert atm_model.species is None
    assert atm_model.vmr is None
    assert atm_model.radius is None


def test_run_atmosphere_calc_p_interp_t_from_pt(tmp_path, reset_jupiter):
    reset = {}
    remove = ['tmodel', 'chemistry', 'radmodel']
    reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', calc_nlayers)
    expected_temperature = np.tile(read_t, calc_nlayers)

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    assert atm_model.species is None
    assert atm_model.vmr is None
    assert atm_model.radius is None


def test_run_atmosphere_read_p_calc_t_from_pt(tmp_path, reset_jupiter):
    reset = {}
    remove = ['nlayers', 'chemistry', 'radmodel']
    reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)
    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', read_nlayers)
    expected_temperature = np.tile(calc_t, read_nlayers)

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    assert atm_model.species is None
    assert atm_model.vmr is None
    assert atm_model.radius is None


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# From atm file only:
def test_run_atmosphere_read_ptqr_from_atm(tmp_path, reset_jupiter):
    reset = {}
    remove = ['nlayers', 'tmodel', 'chemistry', 'radmodel']
    reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', read_nlayers)
    expected_temperature = np.tile(read_t, read_nlayers)
    expected_species = ['H2', 'He', 'H2O', 'CO']
    expected_vmr = np.array([0.85, 0.149, 1.0e-4, 1.0e-4])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)
    np.testing.assert_allclose(atm_model.radius, expected_read_radius)


def test_run_atmosphere_read_ptq_calc_r_from_atm(tmp_path, reset_jupiter):
    reset = {}
    remove = ['nlayers', 'tmodel', 'chemistry']
    reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', read_nlayers)
    expected_temperature = np.tile(read_t, read_nlayers)
    expected_species = ['H2', 'He', 'H2O', 'CO']
    expected_vmr = np.array([0.85, 0.149, 1.0e-4, 1.0e-4])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)
    np.testing.assert_allclose(atm_model.radius, expected_read_ptq_radius)


def test_run_atmosphere_read_p_interp_tqr_from_atm(tmp_path, reset_jupiter):
    reset = {}
    remove = ['tmodel', 'chemistry', 'radmodel']
    reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', calc_nlayers)
    expected_temperature = np.tile(read_t, calc_nlayers)
    expected_species = ['H2', 'He', 'H2O', 'CO']
    expected_vmr = np.array([0.85, 0.149, 1.0e-4, 1.0e-4])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)
    np.testing.assert_allclose(atm_model.radius, expected_interp_tq_radius)


def test_run_atmosphere_calc_p_interp_tq_calc_r_from_atm(
        tmp_path, reset_jupiter,
    ):
    reset = {}
    remove = ['tmodel', 'chemistry']
    reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', calc_nlayers)
    expected_temperature = np.tile(read_t, calc_nlayers)
    expected_species = ['H2', 'He', 'H2O', 'CO']
    expected_vmr = np.array([0.85, 0.149, 1.0e-4, 1.0e-4])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)
    np.testing.assert_allclose(atm_model.radius, expected_interp_tq_calc_radius)


def test_run_atmosphere_take_species_from_atm(tmp_path, reset_jupiter):
    reset = {}
    # Read P and T
    remove = ['nlayers', 'tmodel', 'species']
    reset['atmfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'
    # Calculate VMRs, but take species from input_atmfile
    reset['chemistry'] = 'tea'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    atm_model = pb.run(cfg)

    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', read_nlayers)
    expected_temperature = np.tile(read_t, read_nlayers)
    expected_species = ['H2', 'He', 'H2O', 'CO']
    expected_vmr = np.array([
        8.58283368e-01, 1.40875555e-01, 3.45814052e-04, 4.95262722e-04,
    ])
    output_vmr = atm_model.vmr[0]

    np.testing.assert_allclose(atm_model.press, expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atm_model.temp, expected_temperature)
    np.testing.assert_equal(atm_model.species, expected_species)
    np.testing.assert_allclose(output_vmr, expected_vmr)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Trigger errors:
@pytest.mark.parametrize('atm_input', ('none', 'pt'))
def test_run_atmosphere_calc_pt_missing_q(capfd, tmp_path, atm_input, reset_jupiter):
    reset = {}
    remove = ['chemistry']
    if atm_input == 'none':
        remove += ['atmfile']
    elif atm_input == 'pt':
        reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    error = re.escape(
        'Cannot compute hydrostatic-equilibrium radius profile.\n'
        'radius model needs to know the composition'
    )
    with pytest.raises(ValueError, match=error):
        pb.run(cfg)


@pytest.mark.parametrize('press', ('read', 'calc'))
@pytest.mark.parametrize('temp', ('read', 'calc'))
def test_run_atmosphere_missing_q_from_pt(capfd, tmp_path, press, temp):
    reset = {}
    remove = ['chemistry']
    if press == 'read':
        remove += ['nlayers']
    if temp == 'read':
        remove += ['tmodel']
    reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    error = re.escape(
        'Cannot compute hydrostatic-equilibrium radius profile.\n'
        'radius model needs to know the composition'
    )
    with pytest.raises(ValueError, match=error):
        pb.run(cfg)


def test_run_atmosphere_missing_p(tmp_path):
    remove = ['atmfile', 'chemistry', 'nlayers']
    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        remove=remove,
    )

    error = re.escape(
        'Cannot compute pressure profile, either set {ptop, pbottom, '
        'nlayers} parameters, or provide an input PT profile (ptfile) '
        'or atmospheric file (atmfile)'
    )
    with pytest.raises(ValueError, match=error):
        pb.run(cfg)


def test_run_atmosphere_calc_p_missing_t(capfd, tmp_path):
    remove = ['atmfile', 'chemistry', 'tmodel']
    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        remove=remove,
    )

    error = re.escape(
        'Cannot compute temperature profile, either set a temperature model '
        '(tmodelname) and parameters (tpars), or provide an input PT '
        'profile (ptfile) or atmospheric file (atmfile)'
    )
    with pytest.raises(ValueError, match=error):
        pb.run(cfg)


@pytest.mark.parametrize('atm_input', ('none', 'pt'))
def test_run_atmosphere_calc_q_missing_species(capfd, tmp_path, atm_input, reset_jupiter):
    reset = {}
    remove = ['species']
    if atm_input == 'none':
        remove += ['atmfile']
    elif atm_input == 'pt':
        reset['ptfile'] = 'inputs/jupiter_isothermal_uniform_vmr.atm'

    cfg = make_config(
        tmp_path,
        'configs/atmosphere_jupiter_calc.cfg',
        reset=reset,
        remove=remove,
    )
    error = re.escape(
        'Cannot compute VMRs. Undefined atmospheric species list (species)'
    )
    with pytest.raises(ValueError, match=error):
        pb.run(cfg)


