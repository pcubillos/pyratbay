import os
import sys
import pytest

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay as pb

os.chdir(ROOT+'tests')


def test_hitran(capfd):
    pb.run('tli_hitran_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file: 'tli_hitran_test.cfg'",
        "Reading input database files:",
        "- Mock_HITRAN_H2O_1.00-1.01um.par",
        "Initial TLI wavelength (um):   1.000  (10000.000 cm-1)",
        "Final   TLI wavelength (um):   1.010  ( 9900.990 cm-1)",
        "There are 1 input database file(s).",
        "Database (1/1): 'HITRAN H2O' (H2O molecule)",
        "Number of temperatures: 294",
        "Number of isotopes: 6",
        "Process HITRAN H2O database between records 0 and 887.",
        "[672 148  62   6]",
        "Writing 888 transition lines.",
        "/HITRAN_H2O_1.00-1.01um_test.tli'.",
        ]
    for cap in caps:
        assert cap in captured.out


def test_two_files_two_databases(capfd):
    pb.run('tli_hitran_two_files_two_db_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file:",
        "'tli_hitran_two_files_two_db_test.cfg'",
        "Reading input database files:",
        "- 01_hit12.par",
        "- 02_4000-4500_HITEMP2010.par",
        "There are 2 input database file(s).",
        "Initial TLI wavelength (um):   2.000  ( 5000.000 cm-1)",
        "Final   TLI wavelength (um):   3.000  ( 3333.333 cm-1)",
        "Database (1/2): 'HITRAN H2O' (H2O molecule)",
        "Number of temperatures: 294",
        "Number of isotopes: 6",
        "Database (2/2): 'HITRAN CO2' (CO2 molecule)",
        "Number of isotopes: 10",
        "Cumulative number of isotopes per database: [0, 6, 16]",
        "Process HITRAN H2O database between records 42,645 and 66,776.",
        "Process HITRAN CO2 database between records 0 and 138,257.",
        "[12296  5132  3973  2568   147 99541  9119 22521   220  6854     2]",
        "Writing 162,373 transition lines.",
        "/HITRAN_H2O_CO2_2.0-3.0um_test.tli'.",
        ]
    for cap in caps:
        assert cap in captured.out


def test_two_files_one_database(capfd):
    pb.run('tli_hitran_two_files_one_db_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file:",
        "'tli_hitran_two_files_one_db_test.cfg'",
        "Reading input database files:",
        "- 02_3750-4000_HITEMP2010.par",
        "- 02_4000-4500_HITEMP2010.par",
        "There are 2 input database file(s).",
        "Initial TLI wavelength (um):   2.000  ( 5000.000 cm-1)",
        "Final   TLI wavelength (um):   3.000  ( 3333.333 cm-1)",
        "Database (1/1): 'HITRAN CO2' (CO2 molecule)",
        "Number of temperatures: 294",
        "Number of isotopes: 10",
        "Cumulative number of isotopes per database: [0, 10]",
        "Process HITRAN CO2 database between records 0 and 213,768.",
        "Process HITRAN CO2 database between records 0 and 138,257.",
        "[267537  36165  39420   1090   7500    314]",
        "Writing 352,026 transition lines.",
        "/HITRAN_CO2_2.0-3.0um_test.tli'.",
        ]
    for cap in caps:
        assert cap in captured.out


def test_exomol(capfd):
    pb.run('tli_exomol_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file: 'tli_exomol_test.cfg'",
        "Reading input database files:",
        "- 14N-1H3__MockBYTe__04999-05000.trans",
        "- 15N-1H3__MockBYTe-15__04999-05000.trans",
        "There are 2 input database file(s).",
        "Initial TLI wavelength (um):   2.000  ( 5000.000 cm-1)",
        "Final   TLI wavelength (um):   2.000  ( 4999.950 cm-1)",
        "Database (1/1): 'Exomol NH3' (NH3 molecule)",
        "Number of temperatures: 2000",
        "Number of isotopes: 2",
        "Process Exomol NH3 database between records 0 and 499.",
        "[500 500]",
        "Writing 1,000 transition lines.",
        "/ExoMol_NH3_2.0-2.00002um_test.tli'.",
        ]
    for cap in caps:
        assert cap in captured.out

