# Copyright (c) 2016-2020 Patricio Cubillos.
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE).

import os

import numpy as np

import pyratbay as pb
import pyratbay.plots as pp

os.chdir(pb.constants.ROOT+'tests')


def posterior_pt():
    pyrat = pb.run('configs/mcmc_transmission_test.cfg', init=True)
    with np.load('MCMC_transmission_test.npz') as d:
        Z  = d['Z']
    itemp = np.arange(np.sum(pyrat.ret.pstep[pyrat.ret.itemp] > 0))
    tpars = pyrat.ret.params[pyrat.ret.itemp]

    posterior = Z[:,itemp]
    tmodel = pyrat.atm.tmodel
    targs  = pyrat.atm.targs
    pressure = pyrat.atm.press
    pp.pt(posterior, pyrat.atm.tmodel, pyrat.atm.targs, tpars, itemp, pressure)

