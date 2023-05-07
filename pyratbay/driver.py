# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'run',
]

import os

import numpy as np
import mc3

from . import atmosphere as pa
from . import constants as pc
from . import io as io
from . import opacity as po
from . import plots as pp
from . import spectrum as ps
from . import tools as pt
from . import Pyrat, Atmosphere


def run(cfile, run_step=None, with_log=True):
    """
    Pyrat Bay initialization driver.

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.
    run_step: String
        DEPRECATED
    with_log: Bool
        Flag to save screen outputs to file (True) or not (False)
        (e.g., to prevent overwritting log of a previous run).
    """
    inputs, log = pt.parse(cfile, with_log)
    runmode = inputs.runmode

    # TBD: deprecate run_step
    if run_step == 'dry':
        return inputs

    # Call lineread:
    if runmode == 'tli':
        if inputs.tlifile is None:
            log.error('Undefined TLI file (tlifile)')
        po.make_tli(
            inputs.dblist, inputs.pflist, inputs.dbtype,
            inputs.tlifile[0], inputs.wllow, inputs.wlhigh,
            inputs.wlunits, log,
        )
        return


    # Initialize atmosphere:
    if runmode == 'atmosphere':
        return Atmosphere(inputs, log)


    # Initialize pyrat object:
    pyrat = Pyrat(inputs, log)

    # Calculate extinction-coefficient file if requested:
    if runmode == 'opacity':
        pyrat.compute_opacity()
        return pyrat

    # Stop and return if requested:
    if run_step == 'init':
        return pyrat

    # Compute spectrum and return pyrat object if requested:
    if runmode == "spectrum":
        pyrat.run()
        return pyrat

    if runmode == 'radeq':
        pyrat.radiative_equilibrium()
        return pyrat

    ret = pyrat.ret
    if runmode == 'mcmc' and ret.mcmcfile is None:
        log.error('Undefined MCMC file (mcmcfile)')
    #if pyrat.inputs.resume: # Bypass writting all of the initialization log:
    #    pyrat = Pyrat(args, log=None)
    #    pyrat.log = log
    #else:
    #    pyrat = Pyrat(args, log=log)

    # Mute logging in pyrat object, but not in mc3:
    pyrat.log = mc3.utils.Log(verb=-1, width=80)
    pyrat.spec.specfile = None  # Avoid writing spectrum file during MCMC
    retmodel = False  # Return only the band-integrated spectrum
    # Basename of the output files (no extension):
    outfile = os.path.splitext(ret.mcmcfile)[0]

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

    atm = pyrat.atm
    header = "# MCMC best-fitting atmospheric model.\n\n"
    bestatm = f"{outfile}_bestfit_atmosphere.atm"
    io.write_atm(
        bestatm, atm.press, atm.temp, atm.species,
        atm.vmr, radius=atm.radius,
        punits=atm.punits, runits='km', header=header)

    pyrat.plot_spectrum(spec='best', filename=f'{outfile}_bestfit_spectrum.png')

    # Temperature profiles:
    if atm.temp_model is not None:
        tparams = atm.tpars
        tparams[ret.map_pars['temp']] = bestp[ret.itemp]
        ret.temp_best = atm.temp_model(tparams)

        nsamples, nfree = np.shape(posterior)
        t_posterior = np.tile(tparams, (nsamples,1))
        # Map temperature free parameters from posterior to tparams:
        ifree = np.where(pyrat.ret.pstep>0)[0]
        for j, imap in zip(ret.itemp, ret.map_pars['temp']):
            if j in ifree:
                ipost = list(ifree).index(j)
                t_posterior[:,imap] = posterior[:,ipost]
        tpost = pa.temperature_posterior(t_posterior, atm.temp_model)
        ret.temp_median = tpost[0]
        ret.temp_post_boundaries = tpost[1:]
        pyrat.plot_temperature(
            filename=f'{outfile}_posterior_temperature_profile.png')

    is_emission = pyrat.od.rt_path in pc.emission_rt
    is_transmission = pyrat.od.rt_path in pc.transmission_rt

    if is_emission:
        contrib = ps.contribution_function(pyrat.od.depth, atm.press, pyrat.od.B)
    elif is_transmission:
        contrib = ps.transmittance(pyrat.od.depth, pyrat.od.ideep)
    bands_idx = [band.idx for band in pyrat.obs.filters]
    bands_response = [band.response for band in pyrat.obs.filters]
    band_cf = ps.band_cf(contrib, bands_response, pyrat.spec.wn, bands_idx)

    path = 'transit' if is_transmission else 'emission'
    pp.contribution(
        band_cf, 1.0/(pyrat.obs.bandwn*pc.um), path, atm.press, atm.radius,
        atm.rtop, filename=f"{outfile}_bestfit_cf.png"
    )

    pyrat.log = log  # Un-mute
    log.msg(
        "\nOutput MCMC posterior results, log, bestfit atmosphere, "
        "and spectrum:"
        f"\n  {outfile}.npz"
        f"\n  {inputs.logfile}"
        f"\n  {bestatm}"
        f"\n  {pyrat.spec.specfile}\n\n"
    )
    return pyrat
