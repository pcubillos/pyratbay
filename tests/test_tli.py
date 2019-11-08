import os
import pyratbay as pb

os.chdir(pb.constants.ROOT+'tests')


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
        "Number of temperatures: 503",
        "Number of isotopes: 9",
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
        "Number of temperatures: 503",
        "Number of isotopes: 9",
        "Database (2/2): 'HITRAN CO2' (CO2 molecule)",
        "Number of isotopes: 13",
        "Cumulative number of isotopes per database: [0, 9, 22]",
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
        "Number of temperatures: 353",
        "Number of isotopes: 13",
        "Cumulative number of isotopes per database: [0, 13]",
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


def test_repack(capfd):
    pb.run('tli_repack_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file: 'tli_repack_test.cfg'",
        "Reading input database files:",
        "- CO2_hitran_2.50-2.52um_repack-0.01_lbl.dat",
        "There are 1 input database file(s).",
        "Initial TLI wavelength (um):   2.500  ( 4000.000 cm-1)",
        "Final   TLI wavelength (um):   2.520  ( 3968.254 cm-1)",
        "Database (1/1): 'repack hitran CO2' (CO2 molecule)",
        "Number of temperatures: 294",
        "Number of isotopes: 13",
        "Process repack hitran CO2 database between records 0 and 7,814.",
        "[7438  129  116  132]",
        "Writing 7,815 transition lines.",
        "/repack_HITEMP_CO2_2.50-2.52um.tli'.",
        ]
    for cap in caps:
        assert cap in captured.out


def test_pands(capfd):
    pb.run('tli_pands_test.cfg')
    captured = capfd.readouterr()
    caps = [
        "Read command-line arguments from configuration file: 'tli_pands_test.cfg'",
        "Reading input database files:",
        "- Mock_h2ofastfix.bin",
        "There are 1 input database file(s).",
        "Initial TLI wavelength (um):   2.500  ( 4000.000 cm-1)",
        "Final   TLI wavelength (um):   2.501  ( 3998.401 cm-1)",
        "Database (1/1): 'Partridge & Schwenke (1997)' (H2O molecule)",
        "Number of temperatures: 5",
        "Number of isotopes: 4",
        "Process P&S H2O database between records 38 and 10,220.",
        "[9625  207  219  132]",
        "Writing 10,183 transition lines.",
        "/pands_H2O_2.500-2.501um_test.tli'.",
        ]
    for cap in caps:
        assert cap in captured.out

