# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

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


@pt.ignore_system_exit
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
    phy = pyrat.phy
    atm = pyrat.atm
    ret = pyrat.ret
    inputs = pyrat.inputs

    if run_step == 'dry':
        return pyrat

    # Call lineread package:
    if pyrat.runmode == 'tli':
        if pyrat.lt.tlifile is None:
            log.error('Undefined TLI file (tlifile).')
        po.make_tli(
            inputs.dblist, inputs.pflist, inputs.dbtype,
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


    if pyrat.runmode == 'atmosphere' or pt.isfile(atm.atmfile) != 1:
        # Compute pressure-temperature profile:
        if pt.isfile(atm.ptfile) == 1:
            log.msg(f"\nReading pressure-temperature file: '{atm.ptfile}'.")
            units, _, pressure, temperature = io.read_atm(atm.ptfile)[2:4]
            pressure *= pt.u(units[0]) # pressure in barye
        else:
            check_pressure(pyrat)
            pressure = pa.pressure(
                atm.ptop, atm.pbottom, atm.nlayers, 'barye', log)
            check_temp(pyrat)
            temperature = pa.temperature(
                atm.tmodelname, pressure, atm.nlayers, log, atm.tpars)

        abundances = None
        species = None
        radius = None
        # Compute volume-mixing-ratio profiles:
        if atm.chemistry is not None or pyrat.runmode != 'atmosphere':
            check_atm(pyrat)
            species = inputs.species
            abundances = pa.abundance(
                pressure, temperature, species, inputs.elements,
                inputs.uniform, atm.atmfile, atm.punits, inputs.xsolar,
                atm.escale, inputs.solar, log)

        # Compute altitude profile:
        if abundances is not None and atm.rmodelname is not None:
            check_altitude(pyrat)

            # Mean molecular mass:
            mean_mass = pa.mean_weight(abundances, species)
            # Altitude profile:
            radius = pyrat.hydro(
                pressure, temperature, mean_mass, phy.gplanet,
                phy.mplanet, atm.refpressure, phy.rplanet)

    # Return atmospheric model if requested:
    if pyrat.runmode == 'atmosphere':
        header = '# Pyrat bay atmospheric model\n'
        # Guess radius units if not defined (by rplanet):
        if radius is not None and atm.runits is None:
            atm.runits = 'rjup' if phy.rplanet > 0.5*pc.rjup else 'rearth'
        if atm.atmfile is not None:
            io.write_atm(
                atm.atmfile, pressure, temperature, species,
                abundances, radius, atm.punits, atm.runits, header=header)
            log.msg(f"Output atmospheric file: '{atm.atmfile}'.")
        return pressure, temperature, abundances, species, radius


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
    pyrat.setup_spectrum()
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
        resume=inputs.resume, savefile=ret.mcmcfile)

    if mc3_out is None:
        log.error("Error in MC3.")

    bestp = mc3_out['bestp']
    CRlo  = mc3_out['CRlo']
    CRhi  = mc3_out['CRhi']
    stdp  = mc3_out['stdp']
    posterior, zchain, zmask = mc3.utils.burn(mc3_out)
    ret.posterior = posterior
    ret.bestp = bestp

    # Best-fitting model:
    pyrat.spec.specfile = f"{outfile}_bestfit_spectrum.dat"
    ret.spec_best, ret.bestbandflux = pyrat.eval(bestp)

    header = "# MCMC best-fitting atmospheric model.\n\n"
    bestatm = f"{outfile}_bestfit_atmosphere.atm"
    io.write_atm(
        bestatm, atm.press, atm.temp, pyrat.mol.name,
        atm.q, radius=atm.radius,
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
            pyrat.od.depth, atm.press, pyrat.od.B)
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
    if pyrat.atm.tmodelname is None:
        pyrat.log.error("Undefined temperature model (tmodel).")
    if pyrat.atm.tpars is None:
        pyrat.log.error("Undefined temperature-model parameters (tpars).")


def check_atm(pyrat):
    """
    Check the input arguments to calculate the atmospheric model.
    """
    atm = pyrat.atm
    log = pyrat.log

    if atm.chemistry is None:
        log.error("Undefined chemistry model (chemistry).")
    if atm.atmfile is None:
        log.error("Undefined atmospheric file (atmfile).")
    if pyrat.inputs.species is None:
        log.error("Undefined atmospheric species list (species).")

    # Uniform-abundances profile:
    if atm.chemistry == 'uniform':
        if pyrat.inputs.uniform is None:
            log.error("Undefined list of uniform volume mixing ratios "
                     f"(uniform) for {atm.chemistry} chemistry model.")
        nuniform = len(pyrat.inputs.uniform)
        nspecies = len(pyrat.inputs.species)
        if nuniform != nspecies:
            pyrat.log.error(f"Number of uniform abundances ({nuniform}) does "
                            f"not match the number of species ({nspecies}).")
        return

    # TEA abundances:
    if atm.chemistry == 'tea':
        if pyrat.inputs.elements is None:
            log.error("Undefined elemental composition list (elements) for "
                     f"{atm.chemistry} chemistry model.")

    pyrat.inputs.solar = pyrat.inputs.get_default(
        'solar', 'Solar-abundance file',
        pc.ROOT+'pyratbay/data/AsplundEtal2009.txt')
    pyrat.inputs.xsolar = pyrat.inputs.get_default(
        'xsolar', 'Solar-metallicity scaling factor',
        1.0, gt=0.0, wflag=True)


def check_altitude(pyrat):
    """Check input arguments to calculate altitude profile."""
    atm = pyrat.atm
    phy = pyrat.phy
    log = pyrat.log

    rad_vars = ['mplanet', 'rplanet', 'gplanet']
    missing = [rvar for rvar in rad_vars
               if getattr(phy,rvar) is None]

    if len(missing) == 2:
        err = f'either {missing[0]} or {missing[1]}'
    elif len(missing) == 3:
        err = 'at least two of mplanet, rplanet, or gplanet'

    if len(missing) > 0:
        log.error('Cannot compute hydrostatic-equilibrium radius profile.  '
                 f'Must\ndefine {err}.')

    if pyrat.atm.refpressure is None:
        log.error('Cannot compute hydrostatic-equilibrium radius profile.  '
            'Undefined reference pressure level (refpressure).')

