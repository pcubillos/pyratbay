# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os

import numpy as np

import pyratbay.atmosphere as pa
import pyratbay.constants as pc


os.chdir(pc.ROOT+'tests')

nlayers = 11
pressure = pa.pressure('1e-7 bar', '100 bar', nlayers)
vmr0 = np.tile(
    [0.71578, 0.13458, 3.4839e-04, 4.2549e-04, 8.0067e-08, 1.9381e-22],
    (nlayers,1),
)


def test_metalequil_base():
    vmr_model = pa.vmr_models.MetalEquil()
    assert vmr_model.name == 'metal_equil'
    assert vmr_model.pnames == ['[M/H]']
    assert vmr_model.texnames == ['[M/H]']
    assert vmr_model.npars == 1
    assert vmr_model.type == 'equil'


def test_scaleequil_base():
    vmr_model = pa.vmr_models.ScaleEquil('[C/H]')
    assert vmr_model.name == 'scale_equil'
    assert vmr_model.pnames == ['[C/H]']
    assert vmr_model.texnames == ['[C/H]']
    assert vmr_model.npars == 1
    assert vmr_model.type == 'equil'
    assert vmr_model.element == 'C'


def test_ratioequil_base():
    vmr_model = pa.vmr_models.RatioEquil('C/O')
    assert vmr_model.name == 'ratio_equil'
    assert vmr_model.pnames == ['C/O']
    assert vmr_model.texnames == ['C/O']
    assert vmr_model.npars == 1
    assert vmr_model.type == 'equil'
    assert vmr_model.elements == ['C', 'O']


def test_isovmr_base():
    vmr_model = pa.vmr_models.IsoVMR('H2O', pressure)
    assert vmr_model.name == 'log_H2O'
    assert vmr_model.pnames == ['log_H2O']
    assert vmr_model.texnames == [r'$\log\ X_{\rm H2O}$']
    assert vmr_model.npars == 1
    assert vmr_model.type == 'free'


def test_isovmr_eval():
    spec = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    bulk = np.array(['H2', 'He'])
    mol_models = [pa.vmr_models.IsoVMR('H2O', pressure)]
    params = [-4.0]
    vmr = pa.vmr_scale(vmr0, spec, mol_models, params, bulk)
    # All H2O abundances set to constant value:
    imol = spec.index(mol_models[0].species)
    np.testing.assert_allclose(vmr[:,imol], np.tile(10**params[0], nlayers))


def test_scalevmr_base():
    vmr_model = pa.vmr_models.SlantVMR('H2O', pressure)

    pnames = ['slope_H2O', 'log_VMR0_H2O', 'log_p0_H2O', 'min_log_H2O', 'max_log_H2O']
    texnames = [
        '$m_{\\rm H2O}$',
        '$\\log\\ X_{\\rm H2O}^{0}$',
        '$\\log\\ p_{\\rm H2O}^{0}$',
        '$\\log\\ X_{\\rm H2O}^{\\rm min}$',
        '$\\log\\ X_{\\rm H2O}^{\\rm max}$',
    ]

    assert vmr_model.name == 'slant_H2O'
    assert vmr_model.pnames == pnames
    assert vmr_model.texnames == texnames
    assert vmr_model.npars == 5
    assert vmr_model.type == 'free'


def test_scalevmr_str(capsys):
    vmr_model = pa.vmr_models.SlantVMR('H2O', pressure)
    print(vmr_model)
    captured = capsys.readouterr()

    expected_captured = """VMR model name: slant_H2O
Number of parameters: 5
Parameters: ['slope_H2O', 'log_VMR0_H2O', 'log_p0_H2O', 'min_log_H2O', 'max_log_H2O']
"""
    assert captured.out == expected_captured



def test_scalevmr_eval():
    spec = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    bulk = np.array(['H2', 'He'])
    imol = spec.index('H2O')
    mol_models = [pa.vmr_models.ScaleVMR('H2O', pressure, vmr0[:,imol])]
    params = [-1.0]
    vmr = pa.vmr_scale(vmr0, spec, mol_models, params, bulk)
    # All H2O abundances scaled by value:
    expected_vmr = vmr0[:,imol]*10**params[0]
    np.testing.assert_allclose(vmr[:,imol], expected_vmr, rtol=1e-7)


def test_slantvmr_base():
    vmr_model = pa.vmr_models.ScaleVMR('H2O', pressure, vmr0)
    assert vmr_model.name == 'scale_H2O'
    assert vmr_model.pnames == ['log_H2O']
    assert vmr_model.texnames == [r'$\log\ X_{\rm H2O}$']
    assert vmr_model.npars == 1
    assert vmr_model.type == 'free'


def test_slantvmr_eval():
    spec = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    bulk = np.array(['H2', 'He'])
    mol_models = pa.vmr_models.SlantVMR('H2O', pressure)
    #       slope, vmr0, log_p0, vmr_min, vmr_max
    params = [0.25, -4.0, 0.0, -np.inf, -4.0]
    vmr = pa.vmr_scale(vmr0, spec, mol_models, params, bulk)

    imol = spec.index(mol_models.species)
    expected_vmr = [
        1.77827941e-06, 2.98538262e-06, 5.01187234e-06, 8.41395142e-06,
        1.41253754e-05, 2.37137371e-05, 3.98107171e-05, 6.68343918e-05,
        1.00000000e-04, 1.00000000e-04, 1.00000000e-04,
    ]
    np.testing.assert_allclose(vmr[:,imol], expected_vmr)


def test_vmr_models_eval_list():
    spec = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    bulk = np.array(['H2', 'He'])
    mol_models = [
        pa.vmr_models.SlantVMR('H2O', pressure),
        pa.vmr_models.IsoVMR('CH4', pressure),
    ]
    params = [
        [0.25, -4.0, 0.0, -np.inf, -4.0],
        -3.3,
    ]
    vmr = pa.vmr_scale(vmr0, spec, mol_models, params, bulk)

    imol = spec.index(mol_models[0].species)
    expected_vmr = [
        1.77827941e-06, 2.98538262e-06, 5.01187234e-06, 8.41395142e-06,
        1.41253754e-05, 2.37137371e-05, 3.98107171e-05, 6.68343918e-05,
        1.00000000e-04, 1.00000000e-04, 1.00000000e-04,
    ]
    np.testing.assert_allclose(vmr[:,imol], expected_vmr)

    imol = spec.index(mol_models[1].species)
    np.testing.assert_allclose(vmr[:,imol], np.tile(10**params[1], nlayers))

