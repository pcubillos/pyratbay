# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pyratbay as pb
import pyratbay.tools as pt
import pyratbay.opacity.partitions as pf

from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


def test_hitran(capfd):
    pb.run('configs/tli_hitran_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file:",
        "'configs/tli_hitran_test.cfg'",
        "Reading input database files:",
        "tests/inputs/Mock_HITRAN_H2O_1.00-1.01um.par",
        "Initial TLI wavelength (um):   1.000 (10000.000 cm-1)",
        "Final TLI wavelength (um):     1.010 ( 9900.990 cm-1)",
        "There are 1 input database file(s).",
        "Database (1/1): 'HITRAN H2O' (H2O molecule)",
        "Number of temperatures: 1201",
        "Number of isotopes: 9",
        "Process HITRAN H2O database between records 0 and 887.",
        "[672  148  62  6]",
        "Writing 888 transition lines between wavenumbers 9901.39 and 9999.96 cm-1 (1.000\n-- 1.010 um)",
        "/HITRAN_H2O_1.00-1.01um_test.tli'.",
    ]
    for cap in caps:
        assert cap in captured.out


def test_two_files_two_databases(capfd):
    pb.run('configs/tli_hitran_two_files_two_db_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file:",
        "'configs/tli_hitran_two_files_two_db_test.cfg'",
        "Reading input database files:",
        "tests/inputs/01_hit12.par",
        "tests/inputs/mock_02_4000-4500_HITEMP2010.par",
        "There are 2 input database file(s).",
        "Initial TLI wavelength (um):   2.400 ( 4166.667 cm-1)",
        "Final TLI wavelength (um):     2.500 ( 4000.000 cm-1)",
        "Database (1/2): 'HITRAN H2O' (H2O molecule)",
        "Number of temperatures: 1201",
        "Number of isotopes: 9",
        "Database (2/2): 'HITRAN CO2' (CO2 molecule)",
        "Number of isotopes: 13",
        "Cumulative number of isotopes per database: [0, 9, 22]",
        "Process HITRAN H2O database between records 58,383 and 60,573.",
        "Process HITRAN CO2 database between records 0 and 847.",
        "[1,018  423  305  442  775  44  17  12]",
        "Writing 3,036 transition lines between wavenumbers 4000.00 and 4166.65 cm-1\n(2.400 -- 2.500 um)",
        "/HITRAN_H2O_CO2_2.4-2.5um_test.tli'.",
    ]
    for cap in caps:
        assert cap in captured.out


def test_two_files_one_database(capfd):
    pb.run('configs/tli_hitran_two_files_one_db_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file:",
        "'configs/tli_hitran_two_files_one_db_test.cfg'",
        "Reading input database files:",
        "tests/inputs/mock_02_3750-4000_HITEMP2010.par",
        "tests/inputs/mock_02_4000-4500_HITEMP2010.par",
        "There are 2 input database file(s).",
        "Initial TLI wavelength (um):   2.490 ( 4016.064 cm-1)",
        "Final TLI wavelength (um):     2.510 ( 3984.064 cm-1)",
        "Database (1/1): 'HITRAN CO2' (CO2 molecule)",
        "Number of temperatures: 1001",
        "Number of isotopes: 13",
        "Cumulative number of isotopes per database: [0, 13]",
        "Process HITRAN CO2 database between records 0 and 830.",
        "Process HITRAN CO2 database between records 0 and 847.",
        "[1,515  98  43  23]",
        "Writing 1,679 transition lines between wavenumbers 3996.76 and 4003.22 cm-1\n(2.498 -- 2.502 um)",
        "/HITRAN_CO2_2.49-2.51um_test.tli'.",
    ]
    for cap in caps:
        assert cap in captured.out


def test_two_files_one_dbtype(capfd):
    pb.run('configs/tli_hitran_two_files_one_dbtype.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file:",
        "'configs/tli_hitran_two_files_one_dbtype.cfg'",
        "Reading input database files:",
        "tests/inputs/01_hit12.par",
        "tests/inputs/mock_02_4000-4500_HITEMP2010.par",
        "There are 2 input database file(s).",
        "Initial TLI wavelength (um):   2.400 ( 4166.667 cm-1)",
        "Final TLI wavelength (um):     2.500 ( 4000.000 cm-1)",
        "Database (1/2): 'HITRAN H2O' (H2O molecule)",
        "Number of temperatures: 1001",
        "Number of isotopes: 9",
        "Database (2/2): 'HITRAN CO2' (CO2 molecule)",
        "Number of isotopes: 13",
        "Cumulative number of isotopes per database: [0, 9, 22]",
        "Process HITRAN H2O database between records 58,383 and 60,573.",
        "Process HITRAN CO2 database between records 0 and 847.",
        "[1,018  423  305  442  775  44  17  12]",
        "Writing 3,036 transition lines between wavenumbers 4000.00 and 4166.65 cm-1\n(2.400 -- 2.500 um)",
        "/HITRAN_H2O_CO2_2.4-2.5um_test.tli'.",
    ]
    for cap in caps:
        assert cap in captured.out


def test_two_files_one_pflist(capfd):
    pb.run('configs/tli_hitran_two_files_one_pflist.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file:",
        "'configs/tli_hitran_two_files_one_pflist.cfg'",
        "Reading input database files:",
        "tests/inputs/01_hit12.par",
        "tests/inputs/mock_02_4000-4500_HITEMP2010.par",
        "There are 2 input database file(s).",
        "Initial TLI wavelength (um):   2.400 ( 4166.667 cm-1)",
        "Final TLI wavelength (um):     2.500 ( 4000.000 cm-1)",
        "Database (1/2): 'HITRAN H2O' (H2O molecule)",
        "Number of temperatures: 1001",
        "Number of isotopes: 9",
        "Database (2/2): 'HITRAN CO2' (CO2 molecule)",
        "Number of isotopes: 13",
        "Cumulative number of isotopes per database: [0, 9, 22]",
        "Process HITRAN H2O database between records 58,383 and 60,573.",
        "Process HITRAN CO2 database between records 0 and 847.",
        "[1,018  423  305  442  775  44  17  12]",
        "Writing 3,036 transition lines between wavenumbers 4000.00 and 4166.65 cm-1\n(2.400 -- 2.500 um)",
        "/HITRAN_H2O_CO2_2.4-2.5um_test.tli'.",
    ]
    for cap in caps:
        assert cap in captured.out


def test_exomol(capfd):
    pb.run('configs/tli_exomol_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file:",
        "'configs/tli_exomol_test.cfg'",
        "Reading input database files:",
        "tests/inputs/14N-1H3__MockBYTe__04999-05000.trans",
        "tests/inputs/15N-1H3__MockBYTe-15__04999-05000.trans",
        "There are 2 input database file(s).",
        "Initial TLI wavelength (um):   2.000 ( 5000.000 cm-1)",
        "Final TLI wavelength (um):     2.000 ( 4999.950 cm-1)",
        "Database (1/1): 'Exomol NH3' (NH3 molecule)",
        "Number of temperatures: 2000",
        "Number of isotopes: 2",
        "Process Exomol NH3 database between records 0 and 499.",
        "[500  500]",
        "Writing 1,000 transition lines between wavenumbers 4999.96 and 5000.00 cm-1\n(2.000 -- 2.000 um)",
        "/ExoMol_NH3_2.0-2.00002um_test.tli'.",
    ]
    for cap in caps:
        assert cap in captured.out


def test_repack(capfd):
    pb.run('configs/tli_repack_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file:",
        "'configs/tli_repack_test.cfg'",
        "Reading input database files:",
        "tests/inputs/CO2_hitran_2.50-2.52um_repack-0.01_lbl.dat",
        "There are 1 input database file(s).",
        "Initial TLI wavelength (um):   2.500 ( 4000.000 cm-1)",
        "Final TLI wavelength (um):     2.520 ( 3968.254 cm-1)",
        "Database (1/1): 'repack hitran CO2' (CO2 molecule)",
        "Number of temperatures: 294",
        "Number of isotopes: 13",
        "Process repack hitran CO2 database between records 0 and 7,814.",
        "[7,438  129  116  132]",
        "Writing 7,815 transition lines between wavenumbers 3970.00 and 4000.00 cm-1\n(2.500 -- 2.519 um).",
        "/repack_HITEMP_CO2_2.50-2.52um.tli'.",
    ]
    for cap in caps:
        assert cap in captured.out


def test_pands(capfd):
    with pt.cd('outputs/'):
        pf.kurucz(
            f'{ROOT}tests/inputs/mock_h2opartfn.dat',
            outfile='default',
        )
    pb.run('configs/tli_pands_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file:",
        "'configs/tli_pands_test.cfg'",
        "Reading input database files:",
        "tests/inputs/mock_h2ofastfix.bin",
        "There are 1 input database file(s).",
        "Initial TLI wavelength (um):   2.500 ( 4000.000 cm-1)",
        "Final TLI wavelength (um):     2.501 ( 3998.401 cm-1)",
        "Database (1/1): 'Partridge & Schwenke (1997)' (H2O molecule)",
        "Number of temperatures: 5",
        "Number of isotopes: 4",
        "Process P&S H2O database between records 38 and 10,220.",
        "[9,625  207  219  132]",
        "Writing 10,183 transition lines between wavenumbers 3998.40 and 4000.00 cm-1\n(2.500 -- 2.501 um)",
        "/pands_H2O_2.500-2.501um_test.tli'.",
    ]
    for cap in caps:
        assert cap in captured.out

