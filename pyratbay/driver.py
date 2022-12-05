# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'run',
]

import os

import numpy as np
import mc3

from . import tools as pt
from . import opacity as po
from . import constants as pc
from . import atmosphere as pa
from . import spectrum as ps
from . import io as io
from . import plots as pp
from .pyrat import Pyrat


#@pt.ignore_system_exit
def run(cfile, run_step='run', no_logfile=False):
    """
    Pyrat Bay initialization driver.

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.
    run_step: String
        If 'dry': only read the configuration file into a Pyrat object
        If 'init': initialize a Pyrat object (but no spectra calculation).
        If 'run': run all the way (default).
    no_logfile: Bool
        If True, enforce not to write outputs to a log file
        (e.g., to prevent overwritting log of a previous run).
    """
    pyrat = Pyrat(cfile, no_logfile=no_logfile)
    log = pyrat.log
    atm = pyrat.atm
    ret = pyrat.ret
    inputs = pyrat.inputs

    if run_step == 'dry':
        return pyrat

    # Call lineread:
    if pyrat.runmode == 'tli':
        if pyrat.lt.tlifile is None:
            log.error('Undefined TLI file (tlifile).')
        po.make_tli(
            inputs.dblist, inputs.pflist, inputs.dbtype,
            pyrat.lt.tlifile[0], pyrat.spec.wllow, pyrat.spec.wlhigh,
            pyrat.spec.wlunits, log,
        )
        return


    # Initialize atmosphere:
    pyrat.set_atmosphere()
    if pyrat.runmode == 'atmosphere':
        return pyrat.atm


    # Check status of extinction-coefficient file if necessary:
    if pyrat.runmode != 'spectrum' and pt.isfile(pyrat.ex.extfile) == -1:
        log.error("Undefined extinction-coefficient file (extfile).")

    # Force to re-calculate extinction-coefficient file if requested:
    if pyrat.runmode == 'opacity' and pt.isfile(pyrat.ex.extfile) == 1:
        for extfile in pyrat.ex.extfile:
            os.remove(extfile)

    if pyrat.runmode == 'mcmc' and ret.mcmcfile is None:
        log.error('Undefined MCMC file (mcmcfile).')

    # Initialize pyrat object:
    pyrat.set_spectrum()
    #if pyrat.inputs.resume: # Bypass writting all of the initialization log:
    #    pyrat = Pyrat(args, log=None)
    #    pyrat.log = log
    #else:
    #    pyrat = Pyrat(args, log=log)

    # Stop and return if requested:
    if run_step == 'init' or pyrat.runmode == 'opacity':
        return pyrat

    # Compute spectrum and return pyrat object if requested:
    if pyrat.runmode == "spectrum":
        pyrat.run()
        return pyrat

    if pyrat.runmode == 'radeq':
        pyrat.radiative_equilibrium()
        return pyrat

    # Mute logging in pyrat object, but not in mc3:
    pyrat.log = mc3.utils.Log(verb=0, width=80)
    pyrat.spec.specfile = None  # Avoid writing spectrum file during MCMC
    retmodel = False  # Return only the band-integrated spectrum
    # Basename of the output files (no path, no extension):
    outfile = os.path.splitext(os.path.basename(ret.mcmcfile))[0]

    # Run MCMC:
    mc3_out = mc3.sample(
        data=pyrat.obs.data, uncert=pyrat.obs.uncert,
        func=pyrat.eval, indparams=[retmodel], params=ret.params,
        pmin=ret.pmin, pmax=ret.pmax, pstep=ret.pstep,
        prior=ret.prior, priorlow=ret.priorlow, priorup=ret.priorup,
        sampler=ret.sampler, nsamples=ret.nsamples,
        nchains=ret.nchains, burnin=ret.burnin, thinning=ret.thinning,
        grtest=True, grbreak=ret.grbreak, grnmin=ret.grnmin,
        log=log, ncpu=pyrat.ncpu,
        plots=True, showbp=True,
        pnames=ret.pnames, texnames=ret.texnames,
        resume=inputs.resume, savefile=ret.mcmcfile,
    )

    if mc3_out is None:
        log.error("Error in MC3.")

    bestp = mc3_out['bestp']
    ret.bestp = bestp
    posterior, zchain, zmask = mc3.utils.burn(mc3_out)
    ret.posterior = posterior

    # Best-fitting model:
    pyrat.spec.specfile = f"{outfile}_bestfit_spectrum.dat"
    ret.spec_best, ret.bestbandflux = pyrat.eval(bestp)

    header = "# MCMC best-fitting atmospheric model.\n\n"
    bestatm = f"{outfile}_bestfit_atmosphere.atm"
    io.write_atm(
        bestatm, atm.press, atm.temp, pyrat.mol.name,
        atm.vmr, radius=atm.radius,
        punits=atm.punits, runits='km', header=header)

    pyrat.plot_spectrum(spec='best', filename=f'{outfile}_bestfit_spectrum.png')

    # Temperature profiles:
    if atm.tmodelname is not None:
        ifree = ret.pstep[ret.itemp] > 0
        itemp = np.arange(np.sum(ifree))

        ret.temp_best = atm.tmodel(bestp[ret.itemp])
        tpost = pa.temperature_posterior(
            posterior[:,itemp], atm.tmodel,
            ret.params[ret.itemp], ifree, atm.press)
        ret.temp_median = tpost[0]
        # 1sigma-low, 1sigma-high, 2sigma-low, 2sigma-high:
        ret.temp_post_boundaries = tpost[1:]
        pyrat.plot_temperature(
            filename=f'{outfile}_posterior_temperature_profile.png')

    is_emission = pyrat.od.rt_path in pc.emission_rt
    is_transmission = pyrat.od.rt_path in pc.transmission_rt

    if is_emission:
        cf = ps.contribution_function(
            pyrat.od.depth, atm.press, pyrat.od.B,
        )
        bcf = ps.band_cf(
            cf, pyrat.obs.bandtrans, pyrat.spec.wn, pyrat.obs.bandidx)
    elif is_transmission:
        transmittance = ps.transmittance(pyrat.od.depth, pyrat.od.ideep)
        bcf = ps.band_cf(
            transmittance, pyrat.obs.bandtrans, pyrat.spec.wn,
            pyrat.obs.bandidx)

    path = 'transit' if is_transmission else 'emission'
    pp.contribution(
        bcf, 1.0/(pyrat.obs.bandwn*pc.um), path, atm.press, atm.radius,
        atm.rtop, filename=f"{outfile}_bestfit_cf.png")

    pyrat.log = log  # Un-mute
    log.msg(
       "\nOutput MCMC posterior results, log, bestfit atmosphere, "
       "and spectrum:"
       f"\n'{outfile}.npz'"
       f"\n'{os.path.basename(inputs.logfile)}'"
       f"\n'{bestatm}'"
       f"\n'{pyrat.spec.specfile}'\n\n")
    return pyrat
