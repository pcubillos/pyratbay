# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import multiprocessing as mp

import numpy as np
import scipy.constants as sc
import scipy.interpolate as si
import scipy.special as ss

from .. import tools as pt
from .. import constants as pc
from .. import spectrum as ps
from .. import atmosphere as pa
from .. import io as io
from . import objects as ob


def check_atmosphere(pyrat):
    """
    Check that user input arguments make sense.
    """
    # Check that input files exist:
    if pyrat.mol.molfile is None:
        pyrat.mol.molfile = pc.ROOT + 'pyratbay/data/molecules.dat'

    with pt.log_error(pyrat.log):
        pt.file_exists('molfile', 'Molecular-data', pyrat.mol.molfile)


def check_spectrum(pyrat):
    """
    Check that user input arguments make sense.
    """
    # Shortcuts:
    log = pyrat.log
    phy = pyrat.phy
    spec = pyrat.spec
    atm = pyrat.atm
    obs = pyrat.obs

    # Check that input files exist:
    if pyrat.mol.molfile is None:
        pyrat.mol.molfile = pc.ROOT + 'pyratbay/data/molecules.dat'

    with pt.log_error(log):
        pt.file_exists('atmfile', 'Atmospheric', atm.atmfile)
        pt.file_exists('tlifile', 'TLI', pyrat.lt.tlifile)
        pt.file_exists('molfile', 'Molecular-data', pyrat.mol.molfile)

    if pyrat.runmode == 'spectrum' and spec.specfile is None:
        log.error('Undefined output spectrum file (specfile).')

    if atm.vmr is None and atm.chemistry is None:
        log.error(
            'Missing atmospheric volume mixing ratios. Need to either read '
            'an input profile or compute one via a chemistry model (chemistry)'
        )

    # Compute the Hill radius for the planet:
    phy.rhill = pa.hill_radius(phy.smaxis, phy.mplanet, phy.mstar)

    # Check that the radius profile exists or can be computed:
    if atm.radius is None and pyrat.runmode != 'opacity':
        log.error(
            'Missing atmospheric radius profile.  Need to either read an '
            'input profile or compute one via the radmodel argument'
        )

    # Check Voigt-profile arguments:
    dmin = pyrat.voigt.dmin
    dmax = pyrat.voigt.dmax
    if dmin is not None and dmax is not None and dmax <= dmin:
        log.error(f'dmax ({dmax:g} cm-1) must be > dmin ({dmin:g} cm-1)')

    lmin = pyrat.voigt.lmin
    lmax = pyrat.voigt.lmax
    if lmin is not None and lmax is not None and lmax <= lmin:
        log.error(f'lmax ({lmax:g} cm-1) must be > lmin ({lmin:g} cm-1)')

    if pyrat.runmode == 'opacity' or pt.isfile(pyrat.ex.extfile) == 0:
        if pyrat.ex.tmin is None:
            log.error(
                'Undefined lower temperature boundary (tmin) for '
                'extinction-coefficient grid',
            )
        if pyrat.ex.tmax is None:
            log.error(
                'Undefined upper temperature boundary (tmax) for '
                'extinction-coefficient grid',
            )
        if pyrat.ex.tstep is None:
            log.error(
                'Undefined temperature sampling step (tstep) for '
                'extinction-coefficient grid',
            )
        if pyrat.lt.tlifile is None:
            log.error(
                'Requested extinction-coefficient table, but there '
                'are no input TLI files',
            )

    if pyrat.runmode == 'mcmc':
        if pyrat.od.rt_path in pc.emission_rt:
            if phy.rplanet is None or phy.rstar is None:
                log.error("Undefined radius ratio (need rplanet and rstar)")
        if obs.data is None:
            log.error("Undefined transit/eclipse data (data)")
        if obs.uncert is None:
            log.error("Undefined data uncertainties (uncert)")
        if obs.filters is None:
            log.error("Undefined transmission filters (filters)")
        if len(pyrat.ret.retflag) == 0:
            log.error('Undefined fitting parameters (retrieval_params)')
        if pyrat.ret.sampler is None:
            log.error(
                'Undefined retrieval algorithm (sampler).  '
                'Select from [snooker]'
            )
        if pyrat.ret.nsamples is None:
            log.error('Undefined number of retrieval samples (nsamples)')
        if pyrat.ret.burnin is None:
            log.error('Undefined number of retrieval burn-in samples (burnin)')
        if pyrat.ret.nchains is None:
            log.error('Undefined number of retrieval parallel chains (nchains)')
        if pyrat.ret.params is None:
            log.error('Undefined retrieval fitting parameters (params)')

    # Initialize cloud models:
    pyrat.cloud = ob.Cloud(
        pyrat.inputs.clouds, pyrat.inputs.cpars, pyrat.inputs.fpatchy, log,
    )

    # Initialize Rayleigh models:
    pyrat.rayleigh = ob.Rayleigh(
        pyrat.inputs.rayleigh, pyrat.inputs.rpars, log,
    )

    # Check alkali arguments:
    if pyrat.alkali.model_names is not None:
        pyrat.alkali.models = [
            pa.alkali.get_model(name, pyrat.alkali.cutoff)
            for name in pyrat.alkali.model_names
        ]

    # Accept ray-path argument:
    if pyrat.runmode in ['spectrum', 'mcmc'] and pyrat.od.rt_path is None:
        log.error(
            "Undefined radiative-transfer observing geometry (rt_path)."
            f"  Select from {pc.rt_paths}"
        )

    # Check system arguments:
    if pyrat.od.rt_path in pc.transmission_rt and phy.rstar is None:
        log.error(
            'Undefined stellar radius (rstar), required for transmission '
            'calculation'
        )

    # Check raygrid:
    if spec.raygrid[0] != 0:
        log.error('First angle in raygrid must be 0.0 (normal to surface)')
    if np.any(spec.raygrid < 0) or np.any(spec.raygrid > 90):
        log.error('raygrid angles must lie between 0 and 90 deg')
    if np.any(np.ediff1d(spec.raygrid) <= 0):
        log.error('raygrid angles must be monotonically increasing')
    # Store raygrid values in radians:
    spec.raygrid *= sc.degree

    # Gauss quadrature integration variables:
    if spec.quadrature is not None:
        qnodes, qweights = ss.p_roots(spec.quadrature)
        spec.qnodes = 0.5*(qnodes + 1.0)
        spec.qweights = 0.5 * qweights

    # Number of datapoints and filters:
    if obs.data is not None:
        obs.ndata = len(obs.data)
    if obs.filters is not None:
        obs.nfilters = len(obs.filters)
    # Number checks:
    if obs.uncert is not None:
        n_uncert = len(obs.uncert)
        if obs.ndata != n_uncert:
            log.error(
                f'Number of data uncertainty values ({n_uncert}) does not '
                f'match the number of data points ({obs.ndata})'
            )

    if obs.filters is not None and obs.ndata > 0 and obs.ndata != obs.nfilters:
        log.error(
            f'Number of filter bands ({obs.nfilters}) does not match the '
            f'number of data points ({obs.ndata})'
        )

    if pyrat.ncpu >= mp.cpu_count():
        n_cpu = mp.cpu_count()
        log.warning(
            f'Number of requested CPUs ({pyrat.ncpu}) is >= than the number '
            f'of available CPUs ({n_cpu}).  Enforced ncpu to {n_cpu-1}'
        )
        pyrat.ncpu = n_cpu - 1
    log.head('Check spectrum done.')


def setup(pyrat):
    """
    Process retrieval variables: bulk, molvars, molpars.
    Process stellar spectrum.
    Process the oberving filter bands.
    """
    # Shortcuts:
    inputs = pyrat.inputs
    phy = pyrat.phy
    ret = pyrat.ret
    atm = pyrat.atm
    obs = pyrat.obs
    log = pyrat.log

    # Setup bulk and variable-abundance species:
    species = pyrat.atm.species = pyrat.mol.name
    # Non-overlapping species:
    if atm.bulk is not None and len(np.setdiff1d(atm.bulk, species)) > 0:
        missing = np.setdiff1d(atm.bulk, species)
        log.error(
            f'These bulk species are not present in the atmosphere: {missing}'
        )

    atm.parse_abundance_parameters(pyrat.inputs.molvars, log)
    atm.molpars = inputs.molpars
    # TBD: Do I want to allow null molpars?
    if len(inputs.molvars) > 0 and len(atm.molpars) > 0:
        #if len(inputs.molpars) == 0:
        #    log.error('molvars is set, but there are no molpars')
        if len(inputs.molvars) != len(inputs.molpars):
            log.error(
                'There should be one molpars for each molvars:\n'
                f'molvars: {atm.mol_pnames}\nmolpars: {atm.molpars}'
            )

    # Overlapping species:
    if atm.bulk is not None:
        free_vmr = species[atm.ifree]
        bulk_free_species = np.intersect1d(atm.bulk, free_vmr)
        if len(bulk_free_species) > 0:
            log.error(
                'These species were marked as both bulk and '
                f'variable-abundance: {bulk_free_species}'
            )

    # Obtain abundance ratios between the bulk species:
    if atm.bulk is not None:
        atm.ibulk = [list(species).index(mol) for mol in atm.bulk]
        atm.bulkratio, atm.invsrat = pa.ratio(pyrat.atm.vmr, atm.ibulk)

    # Read stellar spectrum model: starspec, kurucz, or blackbody
    if phy.starspec is not None:
        starflux, starwn, star_temps = io.read_spectra(phy.starspec)
    elif phy.kurucz is not None:
        if phy.tstar is None:
            log.error(
                'Undefined stellar temperature (tstar), required for '
                'Kurucz model'
            )
        if phy.log_gstar is None:
            log.error(
                'Undefined stellar gravity (log_gstar) needed for Kurucz model'
            )
        starflux, starwn, kuruczt, kuruczg = ps.read_kurucz(
            phy.kurucz, phy.tstar, phy.log_gstar,
        )
        log.msg(
            f'Input stellar params: T={phy.tstar:7.1f} K, log(g)={phy.log_gstar:4.2f}\n'
            f'Best Kurucz match:    T={kuruczt:7.1f} K, log(g)={kuruczg:4.2f}'
        )
    elif phy.marcs is not None:
        pass
    elif phy.phoenix is not None:
        pass
    elif phy.tstar is not None:
        starwn = pyrat.spec.wn
        starflux = ps.bbflux(starwn, phy.tstar)
    else:
        starflux, starwn = None, None

    phy.starflux = starflux
    phy.starwn = starwn

    # Store interpolated stellar spectrum:
    if phy.starflux is not None:
        # 1D spectra
        if np.ndim(phy.starflux) == 1:
            sinterp = si.interp1d(phy.starwn, phy.starflux)
            pyrat.spec.starflux = sinterp(pyrat.spec.wn)
        # 2D spectra
        else:
            sinterp = si.interp1d(phy.starwn, phy.starflux, axis=1)
            starflux = sinterp(pyrat.spec.wn)
            pyrat.spec.flux_interp = si.interp1d(star_temps, starflux, axis=0)
            pyrat.spec.starflux = pyrat.spec.flux_interp(phy.tstar)

    is_emission = pyrat.od.rt_path in pc.emission_rt
    if pyrat.runmode=='mcmc' and is_emission and starflux is None:
        log.error(
            'Undefined stellar flux model.  Set starspec, kurucz, or '
            'tstar (for a blackbody spectrum)'
        )

    # Set observational variables (for given filters and other parameters):
    if obs.filters is not None:
        pyrat.set_filters()

    # Temperature models and arguments:
    if atm.tmodelname not in pc.tmodels:
        ntemp = 0
        tpnames, ttexnames = [], []
    else:
        atm.tmodel = pa.tmodels.get_model(
            atm.tmodelname,
            pressure=pyrat.atm.press,
            nlayers=pyrat.atm.nlayers,
        )
        ntemp = atm.tmodel.npars
        tpnames = atm.tmodel.pnames
        ttexnames = atm.tmodel.texnames

    setup_retrieval_parameters_retflag(pyrat)

    # TeX unit conversions for masses and radii:
    utex = {
        'mjup': r'$M_{\rm Jup}$',
        'mearth': r'$M_{\oplus}$',
        'kg': 'kg',
        'gram': 'g',
        'rjup': r'$R_{\rm Jup}$',
        'rearth': r'$R_{\oplus}$',
        'km': 'km',
        'm': 'm',
        'cm': 'cm',
    }

    # Retrieval variables:
    if ret.params is not None and len(ret.params) != ret.nparams:
        nparams = len(ret.params)
        log.error(
            f'The number of input fitting parameters (params, {nparams}) does '
            f'not match the number of required parameters ({ret.nparams})'
        )
    if pyrat.runmode == 'mcmc':
        if ret.pstep is None:
            log.error('Missing pstep argument, required for MCMC runs')

    # Check for non-retrieval model/parameters:
    if (pyrat.rayleigh.models != []
        and (pyrat.runmode != 'mcmc' or 'ray' not in ret.retflag)
        and pyrat.rayleigh.pars is None):
        log.error('Rayleigh parameters (rpars) have not been specified')
    if (pyrat.cloud.models != []
        and (pyrat.runmode != 'mcmc' or 'cloud' not in ret.retflag)
        and pyrat.cloud.pars is None):
        log.error('Cloud parameters (cpars) have not been specified')


def setup_retrieval_parameters_retflag(pyrat):
    """
    Check and setup retrieval parameters via retflag argument.
    To be deprecated, use retrieval_parameters instead.
    """
    log = pyrat.log
    atm = pyrat.atm
    ret = pyrat.ret
    retflag = pyrat.ret.retflag

    if 'temp' in retflag and pyrat.atm.tmodelname is None:
        log.error('Requested temp in retflag, but there is no tmodel')
    if 'mol' in retflag:
        if pyrat.inputs.molvars == []:
            log.error("Requested mol in retflag, but there is no 'molvars'")
        # TBD: This will break for pure eq-chem runs
        if pyrat.atm.bulk is None:
            log.error('Requested mol in retflag, but there are no bulk species')
    if 'ray' in retflag and pyrat.rayleigh.models == []:
        log.error('Requested ray in retflag, but there are no rayleigh models')
    if 'cloud' in retflag and pyrat.cloud.models == []:
        log.error('Requested cloud in retflag, but there are no cloud models')

    # TeX unit conversions for masses and radii:
    utex = {
        'mjup': r'$M_{\rm Jup}$',
        'mearth': r'$M_{\oplus}$',
        'kg': 'kg',
        'gram': 'g',
        'rjup': r'$R_{\rm Jup}$',
        'rearth': r'$R_{\oplus}$',
        'km': 'km',
        'm': 'm',
        'cm': 'cm',
    }

    # Indices to parse the array of fitting parameters:
    ret.nparams = 0
    ret.pnames = []
    ret.texnames = []
    if 'temp' in retflag:
        ntemp = pyrat.atm.tmodel.npars
        ret.itemp = np.arange(ret.nparams, ret.nparams + ntemp)
        ret.pnames += pyrat.atm.tmodel.pnames
        ret.texnames += pyrat.atm.tmodel.texnames
        ret.nparams += ntemp
    if 'rad' in retflag:
        ret.irad = np.arange(ret.nparams, ret.nparams + 1)
        ret.pnames += ['R_planet']
        ret.texnames += [fr'$R_{{\rm p}}$ ({utex[pyrat.atm.runits]})']
        ret.nparams += 1
    if 'press' in retflag:
        ret.ipress = np.arange(ret.nparams, ret.nparams + 1)
        ret.pnames += ['log(p_ref)']
        ret.texnames += [r'$\log p_{{\rm ref}}$']
        ret.nparams += 1
    if 'mol' in retflag:
        nabund = len(atm.mol_pnames)
        ret.imol = np.arange(ret.nparams, ret.nparams + nabund)
        ret.pnames += atm.mol_pnames
        ret.texnames += atm.mol_texnames
        ret.nparams += nabund
    if 'ray' in retflag:
        nray = pyrat.rayleigh.npars
        ret.iray = np.arange(ret.nparams, ret.nparams + nray)
        ret.pnames += pyrat.rayleigh.pnames
        ret.texnames += pyrat.rayleigh.texnames
        ret.nparams += nray
    if 'cloud' in retflag:
        ncloud = pyrat.cloud.npars
        ret.icloud = np.arange(ret.nparams, ret.nparams + ncloud)
        ret.pnames += pyrat.cloud.pnames
        ret.texnames += pyrat.cloud.texnames
        ret.nparams += ncloud
    if 'patchy' in retflag:
        ret.ipatchy = np.arange(ret.nparams, ret.nparams + 1)
        ret.pnames += ['f_patchy']
        ret.texnames += [r'$\phi_{\rm patchy}$']
        ret.nparams += 1
    if 'mass' in retflag:
        ret.imass = np.arange(ret.nparams, ret.nparams + 1)
        ret.pnames += ['M_planet']
        ret.texnames += [fr'$M_{{\rm p}}$ ({utex[pyrat.phy.mpunits]})']
        ret.nparams += 1
    if 'tstar' in retflag:
        ret.itstar = np.arange(ret.nparams, ret.nparams + 1)
        ret.pnames   += ['T_eff']
        ret.texnames += [r'$T_{\rm eff}$ (K)']
        ret.nparams += 1
    if 'offset' in retflag:
        n_offset = len(ret.inst_offset)
        ret.ioffset = np.arange(ret.nparams, ret.nparams + n_offset)
        ret.pnames   += list(ret.inst_offset)
        ret.texnames += list(ret.inst_offset)
        ret.nparams += n_offset

        band_names = [band.name for band in pyrat.obs.filters]
        offset_indices = []
        for inst in ret.inst_offset:
            flags = [inst in name for name in band_names]
            offset_indices.append(flags)
        pyrat.obs.offset_indices = offset_indices

