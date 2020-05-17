# Copyright (c) 2016-2020 Patricio Cubillos.
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE).

__all__ = ["run"]

import os

import numpy as np
import mc3

from . import tools as pt
from . import lineread as lr
from . import constants as pc
from . import atmosphere as pa
from . import spectrum as ps
from . import io as io
from . import plots as pp
from .pyrat import Pyrat


@pt.ignore_system_exit
def run(cfile, init=False, no_logfile=False):
    """
    Pyrat Bay initialization driver.

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.
    init: Bool
        If True, only initialize a Pyrat object (no spectra calculation).
        This is useful when computing spectra interactively.
    no_logfile: Bool
        If True, enforce not to write outputs to a log file
        (e.g., to prevent overwritting log of a previous run).
    """
    pyrat = Pyrat(cfile, no_logfile=no_logfile)
    log = pyrat.log
    phy = pyrat.phy
    atm = pyrat.atm
    ret = pyrat.ret
    inputs = pyrat.inputs

    # Call lineread package:
    if pyrat.runmode == 'tli':
        if pyrat.lt.tlifile is None:
            log.error('Undefined TLI file (tlifile).')
        lr.makeTLI(inputs.dblist, inputs.pflist, inputs.dbtype,
            pyrat.lt.tlifile[0], pyrat.spec.wllow, pyrat.spec.wlhigh,
            pyrat.spec.wlunits, log)
        return


    # Get gplanet from mplanet and rplanet if necessary:
    mplanet = phy.mplanet is not None
    gplanet = phy.gplanet is not None
    rplanet = phy.rplanet is not None

    # Check planetary surface gravity/mass/radius:
    if mplanet and rplanet and gplanet:
        gplanet = pc.G * phy.mplanet / phy.rplanet**2
        if np.abs(gplanet-phy.gplanet)/phy.gplanet > 0.05:
            log.error(
                "All mplanet, rplanet, and gplanet were provided, but "
               f"values are inconsistent (>5%): g(M,R) = {gplanet:7.1f} "
               f"cm s-2 and gplanet = {phy.gplanet:7.1f} cm s-2.")
    elif not mplanet and rplanet and gplanet:
        phy.mplanet = phy.gplanet * phy.rplanet**2 / pc.G
    elif mplanet and not rplanet and gplanet:
        phy.rplanet = np.sqrt(pc.G * phy.mplanet / phy.gplanet)
    elif mplanet and rplanet and not gplanet:
        phy.gplanet = pc.G * phy.mplanet / phy.rplanet**2


    # Compute pressure-temperature profile:
    if pyrat.runmode in ['pt', 'atmosphere'] or pt.isfile(atm.atmfile) != 1:
        if pt.isfile(atm.ptfile) == 1:
            log.msg(f"\nReading pressure-temperature file: '{atm.ptfile}'.")
            pressure, temperature = io.read_pt(atm.ptfile)
        else:
            check_pressure(pyrat)
            pressure = pa.pressure(atm.ptop, atm.pbottom, atm.nlayers,
                'barye', log)
            check_temp(pyrat)
            temperature = pa.temperature(atm.tmodelname, pressure,
                 atm.nlayers, log, atm.tpars)

    # Return temperature-pressure if requested:
    if pyrat.runmode == 'pt':
        return pressure, temperature


    # Compute atmospheric abundances:
    if pyrat.runmode == 'atmosphere' or pt.isfile(atm.atmfile) != 1:
        check_atm(pyrat)
        abundances = pa.abundance(pressure, temperature, inputs.species,
            inputs.elements, inputs.uniform, atm.atmfile, atm.punits,
            inputs.xsolar, atm.escale, inputs.solar, log)

    # Return atmospheric model if requested:
    if pyrat.runmode == 'atmosphere':
        return pressure, temperature, abundances

    # Check status of extinction-coefficient file if necessary:
    if pyrat.runmode != 'spectrum' and pt.isfile(pyrat.ex.extfile) == -1:
        log.error("Undefined extinction-coefficient file (extfile).")

    # Force to re-calculate extinction-coefficient file if requested:
    if pyrat.runmode == 'opacity' and pt.isfile(pyrat.ex.extfile) == 1:
        for extfile in pyrat.ex.extfile:
            os.remove(extfile)

    if pyrat.runmode == 'mcmc' and pyrat.ret.mcmcfile is None:
        log.error('Undefined MCMC file (mcmcfile).')

    # Initialize pyrat object:
    pyrat.setup_spectrum()
    #if pyrat.inputs.resume: # Bypass writting all of the initialization log:
    #    pyrat = Pyrat(args, log=None)
    #    pyrat.log = log
    #else:
    #    pyrat = Pyrat(args, log=log)

    # Stop and return if requested:
    if init or pyrat.runmode == 'opacity':
        return pyrat

    # Compute spectrum and return pyrat object if requested:
    if pyrat.runmode == "spectrum":
        pyrat.run()
        return pyrat


    muted_log = mc3.utils.Log(verb=0, width=80)
    pyrat.log = muted_log      # Mute logging in PB, but not in MC3
    pyrat.spec.specfile = None  # Avoid writing spectrum file during MCMC
    retmodel = False  # Return only the band-integrated spectrum
    # Basename of the output files (no path, no extension):
    outfile = os.path.splitext(os.path.basename(pyrat.ret.mcmcfile))[0]

    # Run MCMC:
    mc3_out = mc3.sample(data=pyrat.obs.data, uncert=pyrat.obs.uncert,
        func=pyrat.eval, indparams=[retmodel], params=ret.params,
        pmin=ret.pmin, pmax=ret.pmax, pstep=ret.pstep,
        prior=ret.prior, priorlow=ret.priorlow, priorup=ret.priorup,
        sampler=ret.sampler, nsamples=ret.nsamples,
        nchains=ret.nchains, burnin=ret.burnin, thinning=ret.thinning,
        grtest=True, grbreak=ret.grbreak, grnmin=ret.grnmin,
        log=log, ncpu=pyrat.ncpu,
        plots=True, showbp=True,
        pnames=ret.pnames, texnames=ret.texnames,
        resume=inputs.resume, savefile=ret.mcmcfile)

    if mc3_out is None:
        log.error("Error in MC3.")

    bestp = mc3_out['bestp']
    CRlo  = mc3_out['CRlo']
    CRhi  = mc3_out['CRhi']
    stdp  = mc3_out['stdp']
    posterior, zchain, zmask = mc3.utils.burn(mc3_out)
    pyrat.ret.posterior = posterior
    pyrat.ret.bestp = bestp

    # Best-fitting model:
    pyrat.spec.specfile = f"{outfile}_bestfit_spectrum.dat"
    pyrat.ret.spec_best, pyrat.ret.bestbandflux = pyrat.eval(bestp)

    header  = "# MCMC best-fitting atmospheric model.\n\n"
    bestatm = f"{outfile}_bestfit_atmosphere.atm"
    io.write_atm(
        bestatm, pyrat.atm.press, pyrat.atm.temp, pyrat.mol.name,
        pyrat.atm.q, radius=pyrat.atm.radius,
        punits=pyrat.atm.punits, runits='km', header=header)

    pyrat.plot_spectrum(spec='best',
        filename=f'{outfile}_bestfit_spectrum.png')

    if pyrat.atm.tmodelname in ['tcea', 'madhu']:
        pyrat.plot_posterior_pt(f'{outfile}_posterior_PT_profile.png')

    if pyrat.od.path == "eclipse":
        cf = ps.contribution_function(
            pyrat.od.depth, pyrat.atm.press, pyrat.od.B)
        bcf = ps.band_cf(
            cf, pyrat.obs.bandtrans, pyrat.spec.wn, pyrat.obs.bandidx)
    elif pyrat.od.path == "transit":
        transmittance = ps.transmittance(pyrat.od.depth, pyrat.od.ideep)
        bcf = ps.band_cf(
            transmittance, pyrat.obs.bandtrans, pyrat.spec.wn,
            pyrat.obs.bandidx)

    pp.contribution(
          bcf, 1.0/(pyrat.obs.bandwn*pc.um), pyrat.od.path,
          pyrat.atm.press, pyrat.atm.radius,
          pyrat.atm.rtop, filename=f"{outfile}_bestfit_cf.png")

    pyrat.log = log  # Un-mute
    log.msg("\nOutput MCMC posterior results, log, bestfit atmosphere, "
        "and spectrum:"
       f"\n'{outfile}.npz'"
       f"\n'{os.path.basename(inputs.logfile)}'"
       f"\n'{bestatm}'"
       f"\n'{pyrat.spec.specfile}'\n\n")
    return pyrat


def check_pressure(pyrat):
    """
    Check the input arguments to calculate the pressure profile.
    """
    if pyrat.atm.nlayers is None:
        pyrat.log.error("Undefined number of atmospheric layers (nlayers).")
    if pyrat.atm.ptop is None:
        pyrat.log.error("Undefined atmospheric top pressure (ptop)")
    if pyrat.atm.pbottom is None:
        pyrat.log.error("Undefined atmospheric bottom pressure (pbottom)")


def check_temp(pyrat):
    """
    Check the input arguments to calculate the temperature profile.
    """
    log = pyrat.log
    atm = pyrat.atm
    if atm.tmodelname is None:
        log.error("Undefined temperature model (tmodel).")
    if atm.tpars is None:
        log.error("Undefined temperature-model parameters (tpars).")


def check_atm(pyrat):
    """
    Check the input arguments to calculate the atmospheric model.
    """
    atm = pyrat.atm
    if atm.atmfile is None:
        pyrat.log.error("Undefined atmospheric file (atmfile).")
    if pyrat.inputs.species is None:
        pyrat.log.error("Undefined atmospheric species list (species).")
    # Uniform-abundances profile:
    if pyrat.inputs.uniform is not None:
        if len(pyrat.inputs.uniform) != len(pyrat.inputs.species):
            pyrat.log.error(
                f"Number of uniform abundances ({len(pyrat.inputs.uniform)}) "
                 "does not match the number of species "
                f"({len(pyrat.inputs.species)}).")
        return
    # TEA abundances:
    if pyrat.inputs.elements is None:
        pyrat.log.error("Undefined atmospheric atomic composition (elements).")
    pyrat.inputs.solar = pyrat.inputs.get_default('solar',
        'Solar-abundance file',
        pc.ROOT+'inputs/AsplundEtal2009.txt')
    pyrat.inputs.atomicfile = pyrat.inputs.get_default('atomicfile',
        'Atomic-composition file', './atomic.tea')
    pyrat.inputs.patm = pyrat.inputs.get_default('patm',
        'Pre-atmospheric file', './preatm.tea')
    pyrat.inputs.xsolar = pyrat.inputs.get_default('xsolar',
        'Solar-metallicity scaling factor', 1.0, gt=0.0, wflag=True)
