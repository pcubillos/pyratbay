# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Loglike',
    'weighted_to_equal',
    'get_multinest_map',
    'multinest_run',
]

import os
import sys

import mc3
import numpy as np
from pymultinest.run import run

from .mpi_tools import (
    get_mpi_rank,
    mpi_barrier,
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
    npars = len(lines) - map_line
    params = np.zeros(npars)
    for i in range(npars):
        params[i] = lines[map_line+i].split()[1]
    return params


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
    os.environ["OMP_NUM_THREADS"] = "1"
    # Synchronize to ensure mkdir call has completed
    mpi_barrier()

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
    output['best_log_post'] = loglike(bestp)
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
