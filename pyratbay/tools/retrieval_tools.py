# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Loglike',
    'weighted_to_equal',
    'get_multinest_map',
    'multinest_run',
    'posterior_post_processing',
]

import os
import sys
import time
import pickle

import mc3
import numpy as np

from ..pyrat import Pyrat
from .. import constants as pc
from .. import plots as pp
from .. import spectrum as ps
from .mpi_tools import get_mpi_rank
from .tools import (
   eta,
   isfile,
)

class Loglike():
    """
    Wrapper to compute the log(likelihood) for a Pyrat object.

    Heavily based on mc3.stats.Loglike class, but this one
    allows to dynamically modify the data and uncertainty values.
    """
    def __init__(self, pyrat):
        self.obs = pyrat.obs
        self.func = pyrat.eval

        self.params = pyrat.ret.params
        self.pstep = pyrat.ret.pstep
        self.ifree = self.pstep>0
        self.ishare = np.where(self.pstep<0)[0]
        # TBD: Throw error when there is no data or retrieval parameters

    def __call__(self, params):
        """
        Evaluate the log(likelihood) for the input set of parameters.

        Parameters
        ----------
        params: 1D float array
            Array of free parameters.
        """
        # Update the free and shared parameters:
        self.params[self.ifree] = params
        for s in self.ishare:
            self.params[s] = self.params[-int(self.pstep[s])-1]

        # Evaluate model (and update data if necessary)
        model = self.func(self.params, retmodel=False)
        data = self.obs.data
        uncert = self.obs.uncert

        log_like = (
            -0.5*np.sum(((data - model) / uncert)**2.0)
            -0.5*np.sum(np.log(2.0*np.pi*uncert**2.0))
        )
        if not np.isfinite(log_like):
            log_like = -1.0e98
        return log_like


def weighted_to_equal(posterior_file, get_weighted=False, min_size=15000):
    """
    Compute an equally-weighted sample from a weighted-probability sample
    read from a Multinest output.

    Parameters
    ----------
    posterior_file: String
        A MultiNest probability-weighted sample output.
    get_weighted: Bool
        If True, also return the weighted sample.
    min_size: Integer
        Set the minimum sample size for the equally weighted posterior.

    Returns
    -------
    equal_posterior: 2D float array
        An equally-weighted posterior sample with dimensions (nsamples, npars).
    weighted_posterior: 2D float array
        The Multinest probabilty-weighted sample with dimensions
        (nsamples, npars).  This is only returned if get_weighted is True.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> posterior = pt.weighted_to_equal('multinest_output.txt')

    >>> # Bet both equal and weighted samples:
    >>> posterior, weighted = pt.weighted_to_equal(
    >>>     'multinest_output.txt',
    >>>     get_weighted=True,
    >>> )
    """
    # MN columns have: sample probability, -2*loglikehood, parameter values
    data = np.loadtxt(posterior_file)
    probability = data[:,0]
    weighted_posterior = data[:,2:]
    nsample = len(probability)

    # Generate PDF from weigths CDF (see Numerical Recipes Sec. 7.3.2)
    # This accomplishes the same as, e.g., dynesty.utils.resample_equal()
    cdf = np.cumsum(probability)
    cdf /= cdf[-1]
    if nsample < min_size:
        nsample = min_size
    rng = np.random.default_rng(seed=None)
    u = sorted(rng.random(nsample))
    indices = np.zeros(nsample, dtype=int)
    i = 0
    i_cdf = 0
    while i < nsample:
        if u[i] < cdf[i_cdf]:
            indices[i] = i_cdf
            i += 1
        else:
            i_cdf += 1
    equal_posterior = weighted_posterior[rng.permutation(indices)]

    if get_weighted:
        return equal_posterior, weighted_posterior
    return equal_posterior


def get_multinest_map(stats_file):
    """
    Get maximum-a-posteriori (MAP) parameters from a MultiNest output file.

    Parameters
    ----------
    stats_file: String
        Path to a Multinest *stats.dat output file.

    Returns
    -------
    params: 1D float array
        The MAP parameter values.
    """
    with open(stats_file, 'r') as f:
        lines = f.readlines()
    map_line = lines.index('MAP Parameters\n') + 2
    nlines = len(lines)

    npars = len(lines) - map_line
    params = []
    for i in range(npars):
        if map_line+i >= nlines:
            break
        if lines[map_line+i].strip() == '':
            break
        index, value = lines[map_line+i].split()
        params.append(value)
    return np.array(params, np.double)


def multinest_run(pyrat, mn_basename):
    """
    A Wrapper of a MultiNest posterior sampling.

    Note
    ----
    For OS X users, it is recommended to set the TMPDIR environment
    variable to "/tmp", e.g., from the command line:
        export TMPDIR=/tmp
    to avoid an MPI error when terminating the execution
    (the call will run to completion in any case)
    https://github.com/open-mpi/ompi/issues/7393#issuecomment-882018321
    """
    from pymultinest.run import run
    os.environ["OMP_NUM_THREADS"] = "1"

    # Shut up for a moment:
    log = pyrat.log
    rank = get_mpi_rank()
    if rank == 0:
        log.msg('Starting Multinest atmospheric retrieval')
    tmp_verb = log.verb
    log.verb = -1

    n_free = np.sum(pyrat.ret.pstep>0)
    prior_transform = mc3.stats.Prior_transform(
        pyrat.ret.prior,
        pyrat.ret.priorlow,
        pyrat.ret.priorup,
        pyrat.ret.pmin,
        pyrat.ret.pmax,
        pyrat.ret.pstep,
    )
    def safe_prior(cube, ndim, nparams):
        try:
            a = np.array([cube[i] for i in range(n_free)])
            b = prior_transform(a)
            for i in range(n_free):
                cube[i] = b[i]
        except Exception as e:
            sys.stderr.write(f'ERROR in prior: {e}\n')
            sys.exit(1)

    loglike = Loglike(pyrat)
    def safe_loglikelihood(cube, ndim, nparams, lnew):
        try:
            a = np.array([cube[i] for i in range(n_free)])
            l = float(loglike(a))
            return l
        except Exception as e:
            sys.stderr.write(f'ERROR in loglikelihood: {e}\n')
            sys.exit(1)


    # The pymultinest call:
    run(
        LogLikelihood=safe_loglikelihood,
        Prior=safe_prior,
        n_dims=n_free,
        importance_nested_sampling=False,
        outputfiles_basename=mn_basename,
        n_live_points=pyrat.ret.nlive,
        resume=pyrat.ret.resume,
        verbose=True,
    )

    if get_mpi_rank() != 0:
        return


    # Post (some plots and stats):
    output = {}
    output['pstep'] = pstep = pyrat.ret.pstep
    output['bestp'] = bestp = pyrat.ret.params
    output['texnames'] = texnames = np.array(pyrat.ret.texnames)
    output['pnames'] = pyrat.ret.pnames
    ifree = np.where(pstep>0)[0]
    ishare = np.where(pstep<0)[0]

    bestp[ifree] = get_multinest_map(f'{mn_basename}stats.dat')
    for s in ishare:
        bestp[s] = bestp[-int(pstep[s])-1]

    posterior, weighted_posterior = weighted_to_equal(
        f'{mn_basename}.txt',
        get_weighted=True,
    )
    output['posterior'] = posterior
    theme = pyrat.ret.theme
    post = mc3.plots.Posterior(
        posterior, pnames=texnames[ifree], theme=theme,
        bestp=bestp[ifree], statistics=pyrat.ret.statistics,
        show_estimates=True,  # TBD: get from cfg?
    )

    # Trace plot:
    savefile = f'{mn_basename}_trace.png'
    mc3.plots.trace(
        weighted_posterior,
        pnames=texnames[ifree],
        color=theme.color,
        savefile=savefile,
    )
    log.msg(savefile, indent=2)
    # Pairwise posteriors plots:
    savefile = f'{mn_basename}_pairwise_posterior.png'
    post.plot(savefile=savefile)
    log.msg(savefile, indent=2)
    # Histogram plots:
    savefile = f'{mn_basename}_marginal_posterior.png'
    post.plot_histogram(savefile=savefile)
    log.msg(savefile, indent=2)


    # Statistics:
    best_spectrum, best_model = pyrat.eval(bestp)
    data = pyrat.obs.data
    uncert = pyrat.obs.uncert
    best_chisq = np.sum((best_model-data)**2 / uncert**2)
    red_chisq = best_chisq / (pyrat.obs.ndata-n_free)
    if pyrat.obs.ndata <= n_free:
        red_chisq = np.nan

    # TBD: need to add log(prior)
    output['best_log_post'] = loglike(bestp[ifree])
    output['best_chisq'] = best_chisq
    output['red_chisq'] = red_chisq
    output['BIC'] = best_chisq + n_free*np.log(pyrat.obs.ndata)
    output['stddev_residuals'] = np.std(best_model-data)

    sample_stats = mc3.stats.calc_sample_statistics(
        post.posterior, bestp, pstep, calc_hpd=True,
    )
    output['medianp'] = sample_stats[0]
    output['meanp'] = sample_stats[1]
    output['stdp'] = sample_stats[2]
    output['median_low_bounds'] = sample_stats[3]
    output['median_high_bounds'] = sample_stats[4]
    output['mode'] = sample_stats[5]
    output['hpd_low_bounds'] = sample_stats[6]
    output['hpd_high_bounds'] = sample_stats[7]

    stats_file = f'{mn_basename}_statistics.txt'
    mc3.stats.summary_stats(post, output, filename=stats_file)

    # Restore verbosity
    log.verb = tmp_verb

    return output


def posterior_post_processing(cfg_file=None, pyrat=None, suffix=''):
    """
    Compute quantities of interest from a retrieval posterior distribution.
    The produced data is stored into a pickle file with root name based
    on the mcmcfile.

    Parameters
    ----------
    cfg_file: String
        A pyratbay config file of a retrieval run (already executed,
        so the parameter posterior files must already exist).
    pyrat: a Pyrat instance
        A pyrat object of an already executed retrieval.
        Used if cfg_file is None.
    """
    if pyrat is None and cfg_file is None:
        raise ValueError(
            "At least one of the input arguments ('cfg_file' or 'pyrat') "
            "must be provided"
        )
    if cfg_file is not None:
        pyrat = Pyrat(cfg_file, log=False, mute=True)

    # Basename of the output files (no extension):
    basename, extension = os.path.splitext(pyrat.ret.mcmcfile)

    if pyrat.ret.sampler == 'multinest':
        if isfile(basename + '.txt') == 0:
            raise ValueError('MultiNest posterior outputs do not exist')
        posterior = weighted_to_equal(basename + '.txt')
    elif pyrat.ret.sampler == 'snooker':
        pass
        # TBD: read from mc3 npz file

    texnames = pyrat.ret.texnames
    theme = pyrat.ret.theme
    post = mc3.plots.Posterior(posterior, texnames, theme=theme)

    pyrat.spec.specfile = None
    nwave = pyrat.spec.nwave
    wl = 1.0 / pyrat.spec.wn / pc.um

    # Unique posterior samples:
    u, uind, uinv = np.unique(
        post.posterior[:,0], return_index=True, return_inverse=True)
    n_unique = len(u)
    print(f'Computing {len(u):d} models for posteriors post-processing')

    # Array of all model parameters (with unique samples)
    u_posterior = np.repeat([pyrat.ret.params], n_unique, axis=0)
    ifree = pyrat.ret.pstep > 0
    u_posterior[:,ifree] = post.posterior[uind]

    is_emission = pyrat.od.rt_path in pc.emission_rt
    is_transmission = pyrat.od.rt_path in pc.transmission_rt

    bands_idx = [band.idx for band in pyrat.obs.filters]
    bands_response = [band.response for band in pyrat.obs.filters]
    nbands = len(bands_response)
    # Evaluate models:
    models = np.zeros((n_unique, nwave))
    band_models = np.zeros((n_unique, nbands))
    temp = np.zeros((n_unique, pyrat.atm.nlayers))
    vmr = np.zeros((n_unique, pyrat.atm.nlayers, pyrat.atm.nmol))
    cf = np.zeros((n_unique, pyrat.atm.nlayers, pyrat.obs.nfilters))
    t0 = time.time()
    for i in range(n_unique):
        models[i], band_models[i] = pyrat.eval(u_posterior[i])
        temp[i] = pyrat.atm.temp
        vmr[i] = pyrat.atm.vmr
        if is_emission:
            contrib = ps.contribution_function(
                pyrat.od.depth, pyrat.atm.press, pyrat.od.B,
            )
        elif is_transmission:
            contrib = ps.transmittance(pyrat.od.depth, pyrat.od.ideep)
        cf[i] = ps.band_cf(contrib, bands_response, pyrat.spec.wn, bands_idx)
        timeleft = eta(time.time()-t0, i+1, n_unique, fmt='.2f')
        if i%3 == 0:
            eta_text = (
                f'{i+1}/{n_unique} samples, '
                f'{100*(i+1)/n_unique:.2f} % done, '
                f'ETA: {timeleft}'
            )
            print(f'{eta_text:80s}', end='\r', flush=True)
    endline = f'{100*(i+1)/n_unique:6.2f} % done'
    print(f'{endline:80s}', flush=True)

    # Spectra posteriors: median -1sigma +1sigma -2sigma +2sigma
    quantiles = np.array([0.5, 0.15865, 0.84135, 0.02275, 0.97725])
    nquantiles = len(quantiles)
    spectrum_posterior = np.zeros((nquantiles,nwave))
    for i in range(nwave):
        msample = models[uinv,i]
        spectrum_posterior[:,i] = np.percentile(msample, 100.0*quantiles)

    temperature_posterior = np.percentile(temp[uinv], 100.0*quantiles, axis=0)
    vmr_posterior = np.percentile(vmr[uinv], 100.0*quantiles, axis=0)
    cf_posterior = cf[uinv]
    cf_median = np.median(cf_posterior, axis=0)

    # Collect spectroscopically active species:
    active_species = []
    for model in pyrat.opacity.models:
        if not hasattr(model, 'species'):
            continue
        if isinstance(model.species, str):
            model_species = [model.species]
        else:
            model_species = list(model.species)
        if model.name == 'H- bound-free/free-free' and 'H-' in pyrat.atm.species:
            model_species.append('H-')
        for spec in model_species:
            if spec not in active_species:
                active_species.append(spec)

    band_models_posterior = np.percentile(
        band_models[uinv,:], 100.0*quantiles, axis=0,
    )

    units = {
        'depth': pyrat.obs.units,
        'flux': 'erg s-1 cm-2 cm',
        'pressure': 'bar',
        'temperature': 'K',
        'wavelength': 'um',
    }

    outputs = {
        'depth_posterior': None,  # placeholder
        'temperature_posterior': temperature_posterior,
        'vmr_posterior': vmr_posterior,
        'band_models_posterior': band_models_posterior,
        'cf_posterior_median': cf_median,
        'pressure': pyrat.atm.press,
        'wl': wl,
        'band_wl': 1.0 / pyrat.obs.bandwn / pc.um,
        'bands_wl': [band.wl for band in pyrat.obs.filters],
        'bands_response': [band.response for band in pyrat.obs.filters],
        'species': pyrat.atm.species,
        'active_species': active_species,
        'starflux': pyrat.spec.starflux,
        'quantiles': quantiles,
        'units': units,
        'path': pyrat.od.rt_path,
        'data': pyrat.obs.data,
        'uncert': pyrat.obs.uncert,
    }

    if is_transmission:
        outputs['depth_posterior'] = spectrum_posterior
    elif is_emission:
        rprs = pyrat.atm.rplanet / pyrat.phy.rstar
        depth = spectrum_posterior / pyrat.spec.starflux * rprs**2.0
        outputs['depth_posterior'] = depth
        outputs['rprs'] = rprs

    post_file = f'{basename}{suffix}_posteriors_info.pickle'
    with open(post_file, 'wb') as handle:
        pickle.dump(outputs, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Now make some plots:
    theme = pyrat.ret.theme
    logxticks = pyrat.inputs.logxticks
    pp.posteriors(post_file, theme=theme, logxticks=logxticks)

    return outputs
