# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import pytest
import pickle

import numpy as np

import pyratbay.io as io
import pyratbay.tools as pt
import pyratbay.constants as pc
import pyratbay.opacity.partitions as pf

os.chdir(pc.ROOT+'tests')


def mock_pf(epf, temp, pf):
    with open(epf, 'w') as f:
        f.write("\n".join(f"{t:7.1f}  {z:.10e}" for t,z in zip(temp,pf)))


def test_get_tips_molname():
    assert pf.get_tips_molname(1) == 'H2O'
    assert pf.get_tips_molname(6) == 'CH4'
    assert pf.get_tips_molname(45) == 'H2'


def test_get_tips_molname_error():
    match = 'TIPS 2021 database does not contain molecule ID: 0'
    with pytest.raises(ValueError, match=match):
        dummy = pf.get_tips_molname(0)


def test_pf_tips():
    expected_temp = np.arange(0, 6001, 5.0)
    expected_temp[0] = 1.0

    with open(f'{pc.ROOT}pyratbay/data/tips_2021.pkl', 'rb') as p:
        expected_pf = pickle.load(p)['H2O']
    with pt.cd('outputs/'):
        pf_data, isotopes, temp = pf.tips('H2O', outfile='default')

    np.testing.assert_equal(temp, expected_temp)
    np.testing.assert_equal(pf_data[6], expected_pf['262'])
    assert isotopes == list(expected_pf.keys())

    pf_read, iso, temp_read = io.read_pf('outputs/PF_tips_H2O.dat')
    np.testing.assert_allclose(pf_read[6,:], expected_pf['262'], rtol=1e-7)
    np.testing.assert_allclose(temp_read, expected_temp, rtol=1e-7)
    assert list(iso) == list(expected_pf.keys())


@pytest.mark.parametrize(
    'outfile',
    [
        'default',
        './PF_tips_HCN.dat',
        f'{pc.ROOT}/tests/outputs/PF_tips_HCN.dat',
    ],
)
def test_pf_tips_outfile(outfile):
    with pt.cd('outputs/'):
        pf_data, isotopes, temp = pf.tips('HCN', outfile=outfile)
    assert 'PF_tips_HCN.dat' in os.listdir('outputs/')
    os.remove('outputs/PF_tips_HCN.dat')


def test_pf_exomol_single():
    epf = f'{pc.ROOT}tests/outputs/14N-1H4__MockBYTe.pf'
    expected_temp = np.arange(1,6)
    expected_pf = np.logspace(0,1,5)
    mock_pf(epf, expected_temp, expected_pf)

    with pt.cd('outputs/'):
        pf_data, isotopes, temp = pf.exomol_pf(epf, outfile='default')
    np.testing.assert_allclose(pf_data[0,:], expected_pf)
    np.testing.assert_allclose(temp, expected_temp)
    assert list(isotopes) == ['41111']

    pf_read, iso, temp_read = io.read_pf('outputs/PF_exomol_NH4.dat')
    np.testing.assert_allclose(pf_read[0,:], expected_pf, rtol=1e-5)
    np.testing.assert_allclose(temp_read, expected_temp, rtol=1e-7)
    assert list(iso) == ['41111']


def test_check_exomol_files_pf():
    files = [
        '1H-12C-14N__Harris.pf',
        '1H-13C-14N__Larner.pf',
    ]
    file_type, molecule, isotopes = pf.check_exomol_files(files)
    assert file_type == 'pf'
    assert molecule == 'HCN'
    assert isotopes == ['124', '134']


def test_check_exomol_files_states():
    files = [
        '1H-12C-14N__Harris.states',
        '1H-13C-14N__Larner.states',
    ]
    file_type, molecule, isotopes = pf.check_exomol_files(files)
    assert file_type == 'states'
    assert molecule == 'HCN'
    assert isotopes == ['124', '134']


def test_check_exomol_files_bz2_states():
    files = [
        '1H-12C-14N__Harris.states.bz2',
        '1H-13C-14N__Larner.states.bz2',
    ]
    file_type, molecule, isotopes = pf.check_exomol_files(files)
    assert file_type == 'states'
    assert molecule == 'HCN'
    assert isotopes == ['124', '134']


def test_check_exomol_files_mix():
    files = [
        '1H-12C-14N__Harris.pf',
        '1H-13C-14N__Larner.states',
    ]
    file_type, molecule, isotopes = pf.check_exomol_files(files)
    assert file_type == ''
    assert molecule == 'HCN'
    assert isotopes == ['124', '134']


def test_check_exomol_files_molecule_mismatch():
    files = [
        '1H-12C-14N__Harris.pf',
        '1H2-16O__POKAZATEL.pf',
    ]
    match = 'All files must correspond to the same molecule'
    with pytest.raises(ValueError, match=match):
        pf.check_exomol_files(files)


def test_pf_exomol_pf_listed_single():
    epf = f'{pc.ROOT}tests/outputs/14N-1H4__MockBYTe.pf'
    expected_temp = np.arange(1,6)
    expected_pf = np.logspace(0,1,5)
    mock_pf(epf, expected_temp, expected_pf)

    with pt.cd('outputs/'):
        pf_data, isotopes, temp = pf.exomol_pf([epf], outfile='default')
    np.testing.assert_allclose(pf_data[0,:], expected_pf)
    np.testing.assert_allclose(temp, expected_temp)
    assert list(isotopes) == ['41111']

    pf_read, iso, temp_read = io.read_pf('outputs/PF_exomol_NH4.dat')
    np.testing.assert_allclose(pf_read[0,:], expected_pf, rtol=1e-5)
    np.testing.assert_allclose(temp_read, expected_temp, rtol=1e-7)
    assert list(iso) == ['41111']


def test_pf_exomol_pf_two():
    epf1 = f'{pc.ROOT}tests/outputs/14N-1H4__MockBYTe.pf'
    epf2 = f'{pc.ROOT}tests/outputs/15N-1H4__MockBYTe.pf'
    expected_temp1 = np.arange(1,6)
    expected_temp2 = np.arange(1,9)
    expected_pf1 = np.logspace(0,1,5)
    expected_pf2 = np.logspace(0,1,8)
    mock_pf(epf1, expected_temp1, expected_pf1)
    mock_pf(epf2, expected_temp2, expected_pf2)

    with pt.cd('outputs/'):
        pf_data, isotopes, temp = pf.exomol_pf([epf1, epf2], outfile='default')
    np.testing.assert_allclose(pf_data[0,:], expected_pf1,     rtol=1e-5)
    np.testing.assert_allclose(pf_data[1,:], expected_pf2[:5], rtol=1e-5)
    np.testing.assert_allclose(temp, expected_temp1)
    assert list(isotopes) == ['41111', '51111']

    pf_read, iso, temp_read = io.read_pf('outputs/PF_exomol_NH4.dat')
    np.testing.assert_allclose(pf_read[0,:], expected_pf1,     rtol=1e-5)
    np.testing.assert_allclose(pf_read[1,:], expected_pf2[:5], rtol=1e-5)
    np.testing.assert_allclose(temp_read, expected_temp1)
    assert list(iso) == ['41111', '51111']


def test_pf_exomol_states():
    exomol_states = [
        'inputs/1H-12C-14N__MockHarris.states',
    ]
    partition, isotopes, temps = pf.exomol_states(
        exomol_states, tmin=100.0, tmax=1500.0, tstep=100.0,
    )
    expected_temps = np.arange(100, 1501, 100)
    expected_pf = np.array([
        6.00000001, 6.00023532, 6.00716014, 6.040937  , 6.12038858,
        6.25410663, 6.44371221, 6.687976  , 6.9850604 , 7.33345747,
        7.7322969 , 8.18138464, 8.68113366, 9.23245386, 9.83662866,
    ])
    np.testing.assert_allclose(partition[0], expected_pf, rtol=1e-5)
    np.testing.assert_allclose(temps, expected_temps)
    assert isotopes == ['124']


def test_pf_kurucz_H2O():
    with pt.cd('outputs/'):
        pf_data, isotopes, temp = pf.kurucz(
            f'{pc.ROOT}tests/inputs/mock_h2opartfn.dat',
            outfile='PF_kurucz_H2O.dat')
    np.testing.assert_equal(temp, np.arange(10, 51, 10.0))
    np.testing.assert_allclose(pf_data,
        np.array([[ 1.328,  3.349,  6.192,  9.417, 12.962],
                  [ 1.33 ,  3.361,  6.217,  9.457, 13.017],
                  [ 1.332,  3.373,  6.24 ,  9.492, 13.066],
                  [ 1.399,  2.945,  5.085,  7.617, 10.476]]))
    np.testing.assert_equal(isotopes,
        np.array(['1H1H16O', '1H1H17O', '1H1H18O', '1H2H16O']))

    pf_read, iso, temp_read = io.read_pf('outputs/PF_kurucz_H2O.dat')
    np.testing.assert_allclose(pf_data,
        np.array([[ 1.328,  3.349,  6.192,  9.417, 12.962],
                  [ 1.33 ,  3.361,  6.217,  9.457, 13.017],
                  [ 1.332,  3.373,  6.24 ,  9.492, 13.066],
                  [ 1.399,  2.945,  5.085,  7.617, 10.476]]))
    np.testing.assert_equal(temp_read, np.arange(10, 51, 10.0))
    np.testing.assert_equal(iso,
        np.array(['1H1H16O', '1H1H17O', '1H1H18O', '1H2H16O']))


def test_pf_kurucz_TiO():
    with pt.cd('outputs/'):
        pf_data, isotopes, temp = pf.kurucz(
            f'{pc.ROOT}tests/inputs/mock_tiopart.dat',
            outfile='default')
    np.testing.assert_allclose(pf_data,
        np.array([[ 28.829,  54.866,  81.572, 110.039, 141.079],
                  [ 28.97 ,  55.151,  82.002, 110.625, 141.834],
                  [ 29.107,  55.425,  82.417, 111.19 , 142.564],
                  [ 29.24 ,  55.692,  82.821, 111.741, 143.273],
                  [ 29.369,  55.95 ,  83.212, 112.273, 143.96 ]]))
    np.testing.assert_equal(temp, np.arange(10, 51, 10.0))
    np.testing.assert_equal(isotopes, np.array(['66','76','86','96','06']))

    pf_read, iso, temp_read = io.read_pf('outputs/PF_kurucz_TiO.dat')
    np.testing.assert_allclose(pf_read,
        np.array([[ 28.829,  54.866,  81.572, 110.039, 141.079],
                  [ 28.97 ,  55.151,  82.002, 110.625, 141.834],
                  [ 29.107,  55.425,  82.417, 111.19 , 142.564],
                  [ 29.24 ,  55.692,  82.821, 111.741, 143.273],
                  [ 29.369,  55.95 ,  83.212, 112.273, 143.96 ]]))
    np.testing.assert_equal(temp_read, np.arange(10, 51, 10.0))
    np.testing.assert_equal(iso, np.array(['66','76','86','96','06']))
