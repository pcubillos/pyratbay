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
        "Reading input database file 'Mock_HITRAN_H2O_1.00-1.01um.par'.",
        "Initial TLI wavelength (um):   1.000  (10000.000 cm-1)",
        "Final   TLI wavelength (um):   1.010  ( 9900.990 cm-1)",
        "['HITRAN H2O']",
        "Number of temperatures: 294",
        "Number of isotopes: 6",
        "Process HITRAN database between records 0 and 887.",
        "[672 148  62   6]",
        "Writing 888 transition lines.",
        "Generated TLI file: './HITRAN_H2O_1.00-1.01um_test.tli'.",
        ]
    for cap in caps:
        assert cap in captured.out


def test_exomol(capfd):
    pb.run('tli_exomol_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Reading input database file '14N-1H3__MockBYTe__04999-05000.trans'.",
        "Reading input database file '15N-1H3__MockBYTe-15__04999-05000.trans'",
        "There are 2 input database file(s).",
        "Initial TLI wavelength (um):   2.000  ( 5000.000 cm-1)",
        "Final   TLI wavelength (um):   2.000  ( 4999.950 cm-1)",
        "['Exomol NH3']",
        "Number of temperatures: 2000",
        "Number of isotopes: 2",
        "Process Exomol database between records 0 and 499.",
        "[500 500]",
        "Writing 1,000 transition lines.",
        "Generated TLI file: './ExoMol_NH3_2.0-2.00002um_test.tli'.",
        ]
    for cap in caps:
        assert cap in captured.out

