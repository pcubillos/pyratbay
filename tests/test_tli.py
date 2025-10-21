# Copyright (c) 2021-2025 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import re
import pytest

import pyratbay as pb
import pyratbay.tools as pt
import pyratbay.opacity.partitions as pf
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


def test_hitran(capfd):
    pb.run('configs/tli_hitran_test.cfg')
    captured = capfd.readouterr()
    cap_lines = [line.strip() for line in captured.out.split('\n')]

    assert cap_lines[7] == "Read command-line arguments from configuration file:"
    assert cap_lines[8] == "'configs/tli_hitran_test.cfg'"
    assert cap_lines[10] == "Reading input database files:"
    assert "tests/inputs/Mock_HITRAN_H2O_1.00-1.01um.par" in cap_lines[11]
    assert cap_lines[12] == "There are 1 input database file(s)."
    assert cap_lines[14] == "Initial wavelength:   1.000 um (10000.000 cm-1)"
    assert cap_lines[15] == "Final wavelength:     1.010 um ( 9900.990 cm-1)"
    assert cap_lines[16] == "There are 1 different database(s)."
    assert cap_lines[19] == "Process HITRAN H2O database between records 0 and 887."
    assert cap_lines[30] == "Database (1/1): 'HITRAN H2O' (H2O molecule)"
    assert cap_lines[31] == "Number of isotopes with line transitions: 4"
    assert cap_lines[33] == "1  '116'        18.011    0.9973173          672"
    assert cap_lines[34] == "2  '118'        20.015    0.0019998          148"
    assert cap_lines[35] == "3  '117'        19.015    0.0003719           62"
    assert cap_lines[36] == "4  '126'        19.017    0.0003107            6"
    assert cap_lines[37] == "Total: 888 line transitions between 1.000 -- 1.010 um"
    assert cap_lines[39] == "Number of temperatures: 1201"
    assert "HITRAN_H2O_1.00-1.01um_test.tli'." in cap_lines[46]


def test_two_files_one_database(capfd):
    pb.run('configs/tli_hitran_two_files_one_db_test.cfg')
    captured = capfd.readouterr()
    cap_lines = [line.strip() for line in captured.out.split('\n')]

    assert cap_lines[8] == "'configs/tli_hitran_two_files_one_db_test.cfg'"
    assert "tests/inputs/mock_02_3750-4000_HITEMP2010.par" in cap_lines[11]
    assert "tests/inputs/mock_02_4000-4500_HITEMP2010.par" in cap_lines[12]
    assert cap_lines[13] == "There are 2 input database file(s)."
    assert cap_lines[15] == "Initial wavelength:   2.490 um ( 4016.064 cm-1)"
    assert cap_lines[16] == "Final wavelength:     2.510 um ( 3984.064 cm-1)"
    assert cap_lines[17] == "There are 1 different database(s)."

    assert cap_lines[20] == "Process HITRAN CO2 database between records 0 and 830."
    assert cap_lines[31] == "Process HITRAN CO2 database between records 0 and 847."
    assert cap_lines[42] == "Database (1/1): 'HITRAN CO2' (CO2 molecule)"
    assert cap_lines[43] == "Number of isotopes with line transitions: 4"
    assert cap_lines[45] == "1  '266'        43.990    0.9842043        1,515"
    assert cap_lines[46] == "2  '366'        44.993    0.0110574           98"
    assert cap_lines[47] == "3  '628'        45.994    0.0039471           43"
    assert cap_lines[48] == "4  '627'        44.994    0.0007340           23"
    assert cap_lines[49] == "Total: 1,679 line transitions between 2.498 -- 2.502 um"
    assert cap_lines[51] == "Number of temperatures: 1001"
    assert "HITRAN_CO2_2.49-2.51um_test.tli" in cap_lines[58]


@pytest.mark.parametrize(
    'config',
    [
        'two_db',
        'one_dbtype',
        'one_pflist',
    ]
)
def test_two_files_two_databases(capfd, config):
    if config == 'two_db':
        cfg = 'configs/tli_hitran_two_files_two_db_test.cfg'
    elif config == 'one_dbtype':
        cfg = 'configs/tli_hitran_two_files_one_dbtype.cfg'
    elif config == 'one_pflist':
        cfg = 'configs/tli_hitran_two_files_one_pflist.cfg'

    pb.run(cfg)
    captured = capfd.readouterr()
    cap_lines = [line.strip() for line in captured.out.split('\n')]

    assert cap_lines[8] == f"{repr(cfg)}"
    assert "tests/inputs/01_hit12.par" in cap_lines[11]
    assert "tests/inputs/mock_02_4000-4500_HITEMP2010.par" in cap_lines[12]
    assert cap_lines[15] == "Initial wavelength:   2.400 um ( 4166.667 cm-1)"
    assert cap_lines[16] == "Final wavelength:     2.500 um ( 4000.000 cm-1)"
    assert cap_lines[17] == "There are 2 different database(s)."

    assert cap_lines[20] == "Process HITRAN H2O database between records 58,383 and 60,573."
    assert cap_lines[31] == "Database (1/2): 'HITRAN H2O' (H2O molecule)"
    assert cap_lines[32] == "Number of isotopes with line transitions: 4"
    assert cap_lines[34] == "1  '116'        18.011    0.9973173        1,018"
    assert cap_lines[35] == "2  '118'        20.015    0.0019998          423"
    assert cap_lines[36] == "3  '117'        19.015    0.0003719          305"
    assert cap_lines[37] == "4  '126'        19.017    0.0003107          442"
    assert cap_lines[38] == "Total: 2,188 line transitions between 2.400 -- 2.500 um"
    assert cap_lines[40] == "Number of temperatures: 1201"

    assert cap_lines[46] == "Process HITRAN CO2 database between records 0 and 847."
    assert cap_lines[57] == "Database (2/2): 'HITRAN CO2' (CO2 molecule)"
    assert cap_lines[58] == "Number of isotopes with line transitions: 4"
    assert cap_lines[60] == "1  '266'        43.990    0.9842043          775"
    assert cap_lines[61] == "2  '366'        44.993    0.0110574           44"
    assert cap_lines[62] == "3  '628'        45.994    0.0039471           17"
    assert cap_lines[63] == "4  '627'        44.994    0.0007340           12"
    assert cap_lines[64] == "Total: 848 line transitions between 2.498 -- 2.500 um"
    assert cap_lines[66] == "Number of temperatures: 1001"
    assert "HITRAN_H2O_CO2_2.4-2.5um_test.tli" in cap_lines[73]


def test_exomol(capfd):
    pb.run('configs/tli_exomol_test.cfg')
    captured = capfd.readouterr()
    cap_lines = [line.strip() for line in captured.out.split('\n')]

    assert cap_lines[8] == "'configs/tli_exomol_test.cfg'"
    assert "tests/inputs/14N-1H3__MockBYTe__04999-05000.trans" in cap_lines[11]
    assert "tests/inputs/15N-1H3__MockBYTe-15__04999-05000.trans" in cap_lines[12]
    assert cap_lines[13] == "There are 2 input database file(s)."
    assert cap_lines[15] == "Initial wavelength:   2.000 um ( 5000.000 cm-1)"
    assert cap_lines[16] == "Final wavelength:     2.000 um ( 4999.950 cm-1)"
    assert cap_lines[17] == "There are 1 different database(s)."
    assert cap_lines[20] == "Process Exomol NH3 database between records 0 and 499."

    assert cap_lines[31] == "Process Exomol NH3 database between records 0 and 499."
    assert cap_lines[42] == "Database (1/1): 'Exomol NH3' (NH3 molecule)"
    assert cap_lines[43] == "Number of isotopes with line transitions: 2"
    assert cap_lines[45] == "1  '4111'       17.027    0.9958716          500"
    assert cap_lines[46] == "2  '5111'       18.024    0.0036613          500"
    assert cap_lines[47] == "Total: 1,000 line transitions between 2.000 -- 2.000 um"
    assert cap_lines[49] == "Number of temperatures: 2000"
    assert "ExoMol_NH3_2.0-2.00002um_test.tli" in cap_lines[54]


def test_repack(capfd):
    pb.run('configs/tli_repack_test.cfg')
    captured = capfd.readouterr()
    cap_lines = [line.strip() for line in captured.out.split('\n')]

    assert cap_lines[8] == "'configs/tli_repack_test.cfg'"
    assert "CO2_hitran_2.50-2.52um_repack-0.01_lbl.dat" in cap_lines[11]
    assert cap_lines[12] == "There are 1 input database file(s)."
    assert cap_lines[14] == "Initial wavelength:   2.500 um ( 4000.000 cm-1)"
    assert cap_lines[15] == "Final wavelength:     2.520 um ( 3968.254 cm-1)"
    assert cap_lines[16] == "There are 1 different database(s)."
    assert cap_lines[19] == "Process repack hitran CO2 database between records 0 and 719."
    assert cap_lines[30] == "Database (1/1): 'repack hitran CO2' (CO2 molecule)"
    assert cap_lines[31] == "Number of isotopes with line transitions: 4"
    assert cap_lines[33] == "1  '266'        43.990    0.9842043          670"
    assert cap_lines[34] == "2  '366'        44.993    0.0110574           30"
    assert cap_lines[35] == "3  '628'        45.994    0.0039471           12"
    assert cap_lines[36] == "4  '627'        44.994    0.0007340            8"
    assert cap_lines[37] == "Total: 720 line transitions between 2.500 -- 2.502 um"
    assert cap_lines[39] == "Number of temperatures: 1001"
    assert "repack_HITEMP_CO2_2.50-2.52um.tli" in cap_lines[46]


def test_exomol_missing_pf_iso():
    error = (
        "No partition functions found for these isotopes of the "
        "NH3 line list: ['4111' '5111']"
    )
    with pytest.raises(ValueError, match=re.escape(error)):
        pb.run('configs/tli_exomol_test_missing_iso.cfg')


def test_repack_bad_lbl_iso(capfd):
    error = "Unrecognized isotope names for CO2 line-list: ['626', '636']"
    with pytest.raises(ValueError) as excinfo:
        pb.run('configs/tli_repack_test_bad_iso.cfg')
    assert error in str(excinfo.value)


@pytest.mark.skip(reason='Deprecated, outdated linelist')
def test_pands(capfd):
    with pt.cd('outputs/'):
        pf.kurucz(
            f'{ROOT}tests/inputs/mock_h2opartfn.dat',
            outfile='default',
        )
    captured = capfd.readouterr()
    pb.run('configs/tli_pands_test.cfg')
    captured = capfd.readouterr()
    cap_lines = [line.strip() for line in captured.out.split('\n')]
    assert cap_lines[8] == "'configs/tli_pands_test.cfg'"
    assert "tests/inputs/mock_h2ofastfix.bin" in cap_lines[11]
    assert cap_lines[12] == "There are 1 input database file(s)."
    assert cap_lines[14] == "Initial wavelength:   2.500 um ( 4000.000 cm-1)"
    assert cap_lines[15] == "Final wavelength:     2.501 um ( 3998.401 cm-1)"
    assert cap_lines[16] == "There are 1 different database(s)."
    assert cap_lines[19] == "Process P&S H2O database between records 38 and 10,220."
    assert cap_lines[30] == "Database (1/1): 'Partridge & Schwenke (1997)' (H2O molecule)"
    assert cap_lines[31] == "Number of isotopes with line transitions: 4"
    assert cap_lines[33] == "1  '116'        18.011    0.9970000        9,625"
    assert cap_lines[34] == "2  '117'        19.015    0.0005080          207"
    assert cap_lines[35] == "3  '118'        20.015    0.0005080          219"
    assert cap_lines[36] == "4  '126'        19.017    0.0019840          132"
    assert cap_lines[37] == "Total: 10,183 line transitions between 2.500 -- 2.501 um"
    assert cap_lines[39] == "Number of temperatures: 5"
    assert "pands_H2O_2.500-2.501um_test.tli" in cap_lines[46]

