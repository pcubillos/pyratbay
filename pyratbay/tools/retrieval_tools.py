# Copyright (c) 2021-2026 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Loglike',
    'weighted_to_equal',
    'posterior_snapshot',
    'get_multinest_map',
    'multinest_run',
    'posterior_post_processing',
]

import datetime
import os
import sys
import time
import pickle

import mc3
import numpy as np

from ..pyrat import Pyrat
from .. import constants as pc
from .. import plots as pp
from .mpi_tools import (
    get_mpi_rank,
    mpi_barrier,
    MPI_Comm,
)
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
        self.ifree = self.pstep > 0
        self.ishare = np.where(self.pstep<0)[0]

        if pyrat.obs.data is None and pyrat.obs.data_hires is None:
            raise ValueError(
                'Attempting to compute a log-likelihood for a model '
                'with no data'
            )
        if np.sum(self.ifree) == 0:
            raise ValueError(
                'Attempting to compute a log-likelihood for a model '
                'with no free parameters'
            )

        self.pnames = np.array(pyrat.ret.texnames)
        self.retrieval_file = pyrat.ret.retrieval_file
        # dt_snapshot hours to seconds
        self._dt_snapshot = pyrat.inputs.dt_retrieval_snapshot * 3600.0
        if self._dt_snapshot > 0:
            self.timer = time.time()

    def get_data(self):
        """
        Concatenate low- and high-resolution data into a single array
        """
        obs = self.obs
        if obs.data is not None and obs.data_hires is not None:
            data = np.concatenate((obs.data, obs.data_hires))
            uncert = np.concatenate((obs.uncert, obs.uncert_hires))
        elif obs.data is not None:
            data = obs.data
            uncert = obs.uncert
        else:
            data = obs.data_hires
            uncert = obs.uncert_hires
        return data, uncert

    def __call__(self, params):
        """
        Evaluate the log(likelihood) for the input set of parameters.

        Parameters
        ----------
        params: 1D float array
            Array of free parameters.
        """
        # TMP plots
        if get_mpi_rank() == 0 and self._dt_snapshot>0:
            now = time.time()
            if now - self.timer >= self._dt_snapshot:
                posterior_snapshot(self.retrieval_file, self.pnames[self.ifree])
                self.timer = now

        # Update the free and shared parameters:
        self.params[self.ifree] = params
        for s in self.ishare:
            self.params[s] = self.params[-int(self.pstep[s])-1]

        # Evaluate model (and update data if necessary)
        model = self.func(self.params, retmodel=False)
        # Concatenate (low-res) data and high-res data arrays
        data, uncert = self.get_data()

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


def posterior_snapshot(retrieval_file, pnames):
    """
    Take a snapshot of a retrieval run, plot the histogram and traces
    of the parameters.
    """
    root, pfile = os.path.split(retrieval_file)
    if not os.path.exists(f'{root}/{pfile}.txt'):
        #print('No posterior file yet')
        return

    equal_posterior, weighted_posterior = weighted_to_equal(
        f'{root}/{pfile}.txt', get_weighted=True,
    )
    with open(f'{root}/{pfile}resume.dat', 'r') as f:
        lines = f.readlines()
    nsamples = int(lines[1].split()[1]) / 1e6
    today = str(datetime.date.today()).replace('-', '_')
    label = f'{nsamples:06.2f}M__{today}'

    if len(np.unique(equal_posterior[:,0])) == 1:
        #print(f'Not enough samples to generate plots ({label}).')
        return

    post = mc3.plots.Posterior(equal_posterior, pnames)
    mc3.plots.trace(
        weighted_posterior,
        pnames=pnames,
        savefile=f'{root}/tmp_trace_{label}.png',
    )
    post.plot_histogram(savefile=f'{root}/tmp_histograms_{label}.png')


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


def multinest_run(pyrat, basename):
    """
    A Wrapper of a MultiNest posterior sampling.

    Parameters
    ----------
    pyrat: Pyrat() object
    basename: String
        Basename for output files. May contain path.
        Should not contain a file extension.

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
        outputfiles_basename=basename,
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

    bestp[ifree] = get_multinest_map(f'{basename}stats.dat')
    for s in ishare:
        bestp[s] = bestp[-int(pstep[s])-1]

    posterior, weighted_posterior = weighted_to_equal(
        f'{basename}.txt',
        get_weighted=True,
    )
    output['posterior'] = posterior
    theme = pyrat.fig.theme
    post = mc3.plots.Posterior(
        posterior, pnames=texnames[ifree], theme=theme,
        bestp=bestp[ifree], statistics=pyrat.ret.statistics,
        show_estimates=True,  # TBD: get from cfg?
    )

    # Trace plot:
    savefile = f'{basename}_posterior_trace.png'
    mc3.plots.trace(
        weighted_posterior,
        pnames=texnames[ifree],
        color=theme.color,
        savefile=savefile,
    )
    log.msg(savefile, indent=2)

    # Statistics:
    best_model = pyrat.eval(bestp, retmodel=False)
    data, uncert = loglike.get_data()

    ndata = len(data)
    best_chisq = np.sum((best_model-data)**2 / uncert**2)
    red_chisq = best_chisq / (ndata-n_free)
    if ndata <= n_free:
        red_chisq = np.nan

    # TBD: need to add log(prior)
    output['best_log_post'] = loglike(bestp[ifree])
    output['best_chisq'] = best_chisq
    output['red_chisq'] = red_chisq
    output['BIC'] = best_chisq + n_free*np.log(ndata)
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

    stats_file = f'{basename}_statistics.txt'
    mc3.stats.summary_stats(post, output, filename=stats_file)

    # Restore verbosity
    log.verb = tmp_verb

    return output


def posterior_post_processing(cfg_file=None, pyrat=None, contributions=None):
    """
    MPI-compatible to compute retrieval posterior quantities of interest
    The produced data is stored into a pickle file with root name based
    on the logfile.

    Parameters
    ----------
    cfg_file: String
        A pyratbay config file of a retrieval run (already executed,
        so the parameter posterior files must already exist).
    pyrat: a Pyrat instance
        A pyrat object of an already executed retrieval.
        Used if cfg_file is None.
    contributions: Bool
        If True, compute and store the posterior spectra for each
        individual absorber.
    """
    if pyrat is None and cfg_file is None:
        raise ValueError(
            "At least one of the input arguments ('cfg_file' or 'pyrat') "
            "must be provided"
        )
    if cfg_file is not None:
        pyrat = Pyrat(cfg_file, log=False, mute=True)
    else:
        pyrat.log.file = None
        pyrat.log.verb = -1

    pyrat.spec.specfile = None
    ifree = pyrat.ret.pstep > 0
    is_eclipse = pyrat.od.rt_path in pc.eclipse_rt
    is_emission = pyrat.od.rt_path in pc.emission_rt
    is_transmission = pyrat.od.rt_path in pc.transmission_rt

    band_wl = pyrat.obs.band_wl
    nwave = pyrat.spec.nwave
    nbands = pyrat.obs.ndata
    nlayers = pyrat.atm.nlayers
    n_tls = pyrat.tls.n_models
    n_offsets = pyrat.obs.depth.n_offsets
    # Individual-absorber contributions leave-one-out (loo) / one-at-a-time (oat)
    per_absorber = contributions in ['loo', 'oat']
    if per_absorber:
        cs_contributions = pyrat.opacity.collect_contributions()
    else:
        cs_contributions = []
    n_absorbers = len(cs_contributions)

    comm = MPI_Comm()
    rank = comm.rank
    size = comm.size

    # Load posteriors
    if rank == 0:
        basename = pyrat.ret.retrieval_file
        if pyrat.ret.sampler == 'multinest':
            post_file = f'{basename}.txt'
            if isfile(post_file) == 0:
                error = f'Posterior file does not exist {repr(post_file)}'
                raise ValueError(error)
            posterior = weighted_to_equal(post_file)
        elif pyrat.ret.sampler == 'snooker':
            mcmc = np.load(basename + '.npz')
            posterior = mc3.utils.burn(mcmc)[0]

        texnames = np.array(pyrat.ret.texnames)
        post = mc3.plots.Posterior(
            posterior, texnames,
            theme=pyrat.fig.theme,
            statistics=pyrat.ret.statistics,
        )
        # All unique parameter samples
        u, uind, uinv = np.unique(
            post.posterior[:,0], return_index=True, return_inverse=True,
        )
        n_unique = len(u)
        u_posterior = np.repeat([pyrat.ret.params], n_unique, axis=0)
        u_posterior[:,ifree] = post.posterior[uind]
        print(f'Computing {len(u):d} models for posterior post-processing')
    else:
        n_unique = None
        u_posterior = None

    n_unique = comm.bcast(n_unique)
    u_posterior = comm.bcast(u_posterior)

    # Allocate shared-memory arrays
    indices = np.array_split(np.arange(n_unique), size)
    models = comm.allocate_shared((n_unique, nwave))
    band_models = comm.allocate_shared((n_unique, nbands))
    temp = comm.allocate_shared((n_unique, nlayers))
    vmr = comm.allocate_shared((n_unique, nlayers, pyrat.atm.nmol))
    cf = comm.allocate_shared((n_unique, nlayers, nbands))
    data = comm.allocate_shared((n_unique, nbands))
    offset = comm.allocate_shared((n_unique, nbands))

    tls_epsilon = comm.allocate_shared((n_unique, n_tls, nwave))
    tls_spectra = comm.allocate_shared((n_unique, n_tls, nwave))
    tls_offset = comm.allocate_shared((n_unique, nbands))
    cs_models = comm.allocate_shared((n_unique, n_absorbers, nwave))

    # Split and evaluate the samples
    indices = np.array_split(np.arange(n_unique), size)
    t0 = time.time()
    for j,i in enumerate(indices[rank]):
        models[i], band_models[i] = pyrat.eval(u_posterior[i])
        temp[i] = pyrat.atm.temp
        vmr[i] = pyrat.atm.vmr
        cf[i] = pyrat.band_contribution()
        data[i] = pyrat.obs.data

        if n_offsets > 0:
            offset[i] = pyrat.obs.inst_offset
        if n_tls > 0:
            tls_epsilon[i] = pyrat.tls.epsilon
            tls_spectra[i] = pyrat.tls.spectrum
            tls_offset[i] = pyrat.tls.band_offset
        for k,absorber in enumerate(cs_contributions):
            # leave-one-out or one-at-a-time contributions
            if contributions == 'loo':
                skip = [absorber]
            elif contributions == 'oat':
                skip = [spec for spec in cs_contributions if spec != absorber]
            cs_models[i,k], _ = pyrat.eval(u_posterior[i], skip=skip)
        if j%5 == 0 and rank==0:
            timeleft = eta(time.time()-t0, size*j+1, n_unique, fmt='.2f')
            eta_text = (
                f'{size*j+1}/{n_unique} samples, '
                f'{100*(size*j+1)/n_unique:.2f} % done, '
                f'ETA: {timeleft}'
            )
            print(f'{eta_text:60s}', flush=True)

    mpi_barrier()
    if rank == 0:
        total_time = f'{(time.time()-t0)/60.0:.2f} min'
        print(f'100.0 % done in {total_time}', flush=True)
    else:
        return

    # Quantiles for all posterior stats: median -1sigma +1sigma -2sigma +2sigma
    quantiles = np.array([0.5, 0.15865, 0.84135, 0.02275, 0.97725])
    nquantiles = len(quantiles)

    spectrum_posterior = np.zeros((nquantiles, nwave))
    cs_contribution_posterior = np.zeros((nquantiles, n_absorbers, nwave))
    tls_posterior = np.zeros((nquantiles, n_tls, nwave))
    tls_spectra_posterior = np.zeros((nquantiles, n_tls, nwave))
    for i in range(nwave):
        sample = models[uinv,i]
        spectrum_posterior[:,i] = np.quantile(sample, quantiles)
        if n_tls > 0:
            sample = tls_epsilon[uinv,:,i]
            tls_posterior[:,:,i] = np.quantile(sample, quantiles, axis=0)
            sample = tls_spectra[uinv,:,i]
            tls_spectra_posterior[:,:,i] = np.quantile(sample, quantiles, axis=0)
        if per_absorber:
            sample = cs_models[uinv,:,i]
            cs_contribution_posterior[:,:,i] = np.quantile(sample, quantiles, axis=0)

    band_models_posterior = np.quantile(band_models[uinv,:], quantiles, axis=0)
    data_posterior = np.quantile(data[uinv], quantiles, axis=0)
    if n_offsets > 0:
        offset_posterior = np.quantile(offset[uinv], quantiles, axis=0)
    if n_tls > 0:
        tls_offset_posterior = np.quantile(tls_offset[uinv], quantiles, axis=0)

    temperature_posterior = np.quantile(temp[uinv], quantiles, axis=0)
    vmr_posterior = np.quantile(vmr[uinv], quantiles, axis=0)
    cf_posterior = cf[uinv]
    cf_median = np.median(cf_posterior, axis=0)

    # Parameter statistics
    stats_1sigma = mc3.stats.calc_sample_statistics(
        post.posterior, pyrat.ret.params, pyrat.ret.pstep, quantile=0.683,
    )
    stats_2sigma = mc3.stats.calc_sample_statistics(
        post.posterior, pyrat.ret.params, pyrat.ret.pstep, quantile=0.9545,
    )
    nfree = np.sum(ifree)
    params_posterior = np.zeros((nquantiles, nfree))
    params_posterior[0] = stats_1sigma[0][ifree]
    params_posterior[1] = stats_1sigma[3][ifree]
    params_posterior[2] = stats_1sigma[4][ifree]
    params_posterior[3] = stats_2sigma[3][ifree]
    params_posterior[4] = stats_2sigma[4][ifree]

    # Collect spectroscopically active species
    active_species = []
    for model in pyrat.opacity.models:
        if not hasattr(model, 'species'):
            continue
        if isinstance(model.species, str):
            model_species = [model.species]
        else:
            model_species = list(model.species)
        # TBD: remove second condition
        if model.name == 'H- continuum' and 'H-' in pyrat.atm.species:
            model_species.append('H-')
        for spec in model_species:
            if spec not in active_species:
                active_species.append(spec)

    if pyrat.od.rt_path == 'f_lambda':
        flux_units = 'W m-2 um-1'
    else:
        flux_units = 'erg s-1 cm-2 cm'
    units = {
        'depth': pyrat.obs.units,
        'flux': flux_units,
        'pressure': 'bar',
        'temperature': 'K',
        'wavelength': 'um',
    }

    outputs = {}
    if is_transmission:
        outputs['depth_posterior'] = spectrum_posterior
    elif is_emission:
        outputs['flux_posterior'] = spectrum_posterior
    elif is_eclipse:
        rprs = pyrat.atm.rplanet / pyrat.atm.rstar
        fplanet = spectrum_posterior * pyrat.spec.starflux / rprs**2.0
        outputs['depth_posterior'] = spectrum_posterior
        outputs['flux_posterior'] = fplanet
        outputs['rprs'] = rprs

    outputs |= {
        'temperature_posterior': temperature_posterior,
        'vmr_posterior': vmr_posterior,
        'band_models_posterior': band_models_posterior,
        'cf_posterior_median': cf_median,
    }
    if per_absorber:
        outputs['absorber_contribution_labels'] = cs_contributions
        outputs['absorber_contribution_posterior'] = cs_contribution_posterior
    if n_tls > 0:
        outputs['tls_posterior'] = tls_posterior
        outputs['tls_spectra_posterior'] = tls_spectra_posterior
        outputs['tls_offset_posterior'] = tls_offset_posterior
        outputs['tls_labels'] = pyrat.tls.models
        outputs['tls_mask'] = pyrat.tls.band_mask
    if n_offsets > 0:
        outputs['offset_posterior'] = offset_posterior
    outputs |= {
        'params_posterior': params_posterior,
        'params_names': np.array(pyrat.ret.pnames)[ifree],
        'params_texnames': texnames[ifree],
        'pressure': pyrat.atm.press,
        'wl': pyrat.spec.wl,
        'band_wl': band_wl,
        'band_half_widths': pyrat.obs.half_widths,
        'species': pyrat.atm.species,
        'active_species': active_species,
        'starflux': pyrat.spec.starflux,
        'quantiles': quantiles,
        'units': units,
        'path': pyrat.od.rt_path,
    }

    if pyrat.obs.data is not None:
        outputs['data_posterior'] = data_posterior
        outputs['data'] = pyrat.obs.depth.data
        outputs['uncert'] = pyrat.obs.uncert
        outputs['band_labels'] = [band.name for band in pyrat.obs.bands]
    if pyrat.obs.data_hires is not None:
        outputs['data_hires'] = pyrat.obs.data_hires
        outputs['uncert_hires'] = pyrat.obs.uncert_hires

    # Figure plotting configs
    outputs['theme'] = pyrat.fig.theme
    outputs['log_wl'] = pyrat.fig.log_wl
    outputs['fig_resolution'] = pyrat.fig.resolution
    outputs['fig_data_color'] = pyrat.fig.data_color

    post_file = f'{basename}_posteriors_info.pickle'
    with open(post_file, 'wb') as handle:
        pickle.dump(outputs, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Now make some plots
    pp.posteriors(post_file)

