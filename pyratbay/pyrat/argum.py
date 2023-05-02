# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import multiprocessing as mp

import numpy as np
import scipy.interpolate as si

from .. import tools as pt
from .. import constants as pc
from .. import spectrum as ps
from .. import atmosphere as pa
from .. import io as io
from . import objects as ob
from .rayleigh import Rayleigh


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

    with pt.log_error(log):
        pt.file_exists('atmfile', 'Atmospheric', atm.atmfile)
        pt.file_exists('tlifile', 'TLI', pyrat.lt.tlifile)
        pt.file_exists('molfile', 'Molecular-data', atm.molfile)

    if pyrat.runmode == 'spectrum' and spec.specfile is None:
        log.error('Undefined output spectrum file (specfile).')

    if atm.vmr is None and atm.chemistry is None:
        log.error(
            'Missing atmospheric volume mixing ratios. Need to either read '
            'an input profile or compute one via a chemistry model (chemistry)'
        )

    # Compute the Hill radius for the planet:
    atm.rhill = pa.hill_radius(atm.smaxis, atm.mplanet, phy.mstar)

    # Radiative-transfer path:
    pyrat.od.rt_path = pyrat.inputs.rt_path

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
            if atm.rplanet is None or phy.rstar is None:
                log.error("Undefined radius ratio (need rplanet and rstar)")
        if obs.data is None:
            log.error("Undefined transit/eclipse data (data)")
        if obs.uncert is None:
            log.error("Undefined data uncertainties (uncert)")
        if obs.nfilters == 0:
            log.error("Undefined transmission filters (filters)")
        if pyrat.inputs.retrieval_params is None and pyrat.ret.retflag is None:
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
    pyrat.rayleigh = Rayleigh(
        pyrat.inputs.rayleigh,
        pyrat.inputs.rpars,
        pyrat.atm.species,
        pyrat.spec.wn,
        log,
        pyrat.cloud,
    )

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
    species = pyrat.atm.species
    # Non-overlapping species:
    if atm.bulk is not None and len(np.setdiff1d(atm.bulk, species)) > 0:
        missing = np.setdiff1d(atm.bulk, species)
        log.error(
            f'These bulk species are not present in the atmosphere: {missing}'
        )

    atm.parse_abundance_parameters(pyrat.inputs.molvars, log)
    atm.molpars = inputs.molpars
    # Allow null molpars for now, check later after retrieval_params
    if atm.mol_npars > 0 and atm.molpars is not None:
        if atm.mol_npars != len(atm.molpars):
            log.error(
                'There should be one molpars value for each molvars variable:\n'
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

    offset_pnames = obs.offset_instruments
    if offset_pnames is None:
        offset_pnames = []

    # TBD: test run with no tmodel
    if atm.tmodelname in pc.tmodels:
        temp_pnames = atm.temp_model.pnames
    else:
        temp_pnames = []

    if pyrat.inputs.retrieval_params is None:
        setup_retrieval_parameters_retflag(pyrat)
        return

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

    # Indices to map free parameters of each model:
    ret.map_pars = map_pars = {
        'temp': [],
        'mol': [],
        'ray': [],
        'cloud': [],
        'offset': [],
    }
    ret.nparams = len(ret.pnames)
    ret.texnames = [None for _ in ret.pnames]
    # Indices to parse the array of fitting parameters:
    itemp = []
    imol = []
    iray = []
    icloud = []
    ioffset = []

    solo_params = [
        'log_p_ref',
        'R_planet',
        'M_planet',
        'f_patchy',
        'T_eff',
    ]
    all_available_params = (
        solo_params +
        temp_pnames +
        atm.mol_pnames +
        pyrat.rayleigh.pnames +
        pyrat.cloud.pnames +
        offset_pnames
    )

    for i,pname in enumerate(ret.pnames):
        if pname == 'log_p_ref':
            ret.ipress = np.array([i])
            ret.texnames[i] = r'$\log p_{{\rm ref}}$'
        elif pname == 'R_planet':
            ret.irad = np.array([i])
            ret.texnames[i] = fr'$R_{{\rm p}}$ ({utex[pyrat.atm.runits]})'
        elif pname == 'M_planet':
            ret.imass = np.array([i])
            ret.texnames[i] = fr'$M_{{\rm p}}$ ({utex[pyrat.phy.mpunits]})'
        elif pname =='f_patchy':
            ret.ipatchy = np.array([i])
            ret.texnames[i] = r'$\phi_{\rm patchy}$'
        elif pname == 'T_eff':
            ret.itstar = np.array([i])
            ret.texnames[i] = r'$T_{\rm eff}$ (K)'

        elif pname in temp_pnames:
            itemp.append(i)
            idx = atm.temp_model.pnames.index(pname)
            map_pars['temp'].append(idx)
            ret.texnames[i] = atm.temp_model.texnames[idx]
        elif pname in atm.mol_pnames:
            imol.append(i)
            idx = atm.mol_pnames.index(pname)
            map_pars['mol'].append(idx)
            ret.texnames[i] = atm.mol_texnames[idx]
        elif pname in pyrat.rayleigh.pnames:
            iray.append(i)
            idx = pyrat.rayleigh.pnames.index(pname)
            map_pars['ray'].append(idx)
            ret.texnames[i] = pyrat.rayleigh.texnames[idx]
        elif pname in pyrat.cloud.pnames:
            icloud.append(i)
            idx = pyrat.cloud.pnames.index(pname)
            map_pars['cloud'].append(idx)
            ret.texnames[i] = pyrat.cloud.texnames[idx]
        elif pname in offset_pnames:
            ioffset.append(i)
            idx = offset_pnames.index(pname)
            map_pars['offset'].append(idx)
            ret.texnames[i] = pname.replace('offset_', r'$\Delta$')
        else:
            log.error(
                f"Invalid retrieval parameter '{pname}'. Possible "
                f"values are:\n{all_available_params}"
            )

    if len(itemp) > 0:
        ret.itemp = itemp
    if len(imol) > 0:
        ret.imol = imol
    if len(icloud) > 0:
        ret.icloud = icloud
    if len(iray) > 0:
        ret.iray = iray
    if len(ioffset) > 0:
        ret.ioffset = ioffset

    # Patch missing parameters if possible:
    patch_temp = (
        atm.tpars is None and
        ret.itemp is not None and
        len(map_pars['temp']) == atm.temp_model.npars
    )
    if patch_temp:
        atm.tpars = np.zeros(atm.temp_model.npars)
        atm.tpars[map_pars['temp']] = ret.params[ret.itemp]

    patch_abundance = (
        atm.molpars is None and
        ret.imol is not None and
        len(map_pars['mol']) == atm.mol_npars
    )
    if patch_abundance:
        atm.molpars = np.zeros(len(atm.mol_pnames))
        atm.molpars[map_pars['mol']] = ret.params[ret.imol]

    patch_rayleigh = (
        pyrat.rayleigh.pars is None and
        ret.iray is not None and
        len(map_pars['ray']) == pyrat.rayleigh.npars
    )
    if patch_rayleigh:
        pyrat.rayleigh.pars = np.zeros(pyrat.rayleigh.npars)
        pyrat.rayleigh.pars[map_pars['ray']] = ret.params[ret.iray]

    patch_cloud = (
        pyrat.cloud.pars is None and
        ret.icloud is not None and
        len(map_pars['cloud']) == pyrat.cloud.npars
    )
    if patch_cloud:
        pyrat.cloud.pars = np.zeros(pyrat.cloud.npars)
        pyrat.cloud.pars[map_pars['cloud']] = ret.params[ret.icloud]

    if atm.tpars is None and atm.temp_model is not None:
        log.error('Not all temperature parameters were defined (tpars)')
    if atm.molpars is None and atm.mol_npars > 0:
        log.error('Not all abundance parameters were defined (molpars)')
    if pyrat.cloud.pars is None and pyrat.cloud.npars > 0:
        log.error('Not all Cloud parameters were defined (cpars)')
    if pyrat.rayleigh.pars is None and pyrat.rayleigh.npars > 0:
        log.error('Not all Rayleigh parameters were defined (rpars)')


def setup_retrieval_parameters_retflag(pyrat):
    """
    Check and setup retrieval parameters via retflag argument.
    To be deprecated, use retrieval_parameters instead.
    """
    log = pyrat.log
    atm = pyrat.atm
    ret = pyrat.ret
    retflag = pyrat.ret.retflag

    if retflag is None:
        return

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

    ret.map_pars = {
        'temp': [],
        'mol': [],
        'ray': [],
        'cloud': [],
        'offset': [],
    }
    # Indices to parse the array of fitting parameters:
    ret.nparams = 0
    ret.pnames = []
    ret.texnames = []
    if 'temp' in retflag:
        ntemp = pyrat.atm.temp_model.npars
        ret.itemp = np.arange(ret.nparams, ret.nparams + ntemp)
        ret.pnames += pyrat.atm.temp_model.pnames
        ret.texnames += pyrat.atm.temp_model.texnames
        ret.nparams += ntemp
        ret.map_pars['temp'] = np.arange(ntemp)
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
        ret.map_pars['mol'] = np.arange(nabund)
    if 'ray' in retflag:
        nray = pyrat.rayleigh.npars
        ret.iray = np.arange(ret.nparams, ret.nparams + nray)
        ret.pnames += pyrat.rayleigh.pnames
        ret.texnames += pyrat.rayleigh.texnames
        ret.nparams += nray
        ret.map_pars['ray'] = np.arange(nray)
    if 'cloud' in retflag:
        ncloud = pyrat.cloud.npars
        ret.icloud = np.arange(ret.nparams, ret.nparams + ncloud)
        ret.pnames += pyrat.cloud.pnames
        ret.texnames += pyrat.cloud.texnames
        ret.nparams += ncloud
        ret.map_pars['cloud'] = np.arange(ncloud)
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
        n_offset = len(ret.offset_instruments)
        ret.ioffset = np.arange(ret.nparams, ret.nparams + n_offset)
        ret.pnames   += list(ret.offset_instruments)
        ret.texnames += list(ret.offset_instruments)
        ret.nparams += n_offset

        band_names = [band.name for band in pyrat.obs.filters]
        offset_indices = []
        for inst in ret.offset_instruments:
            flags = [inst in name for name in band_names]
            offset_indices.append(flags)
        pyrat.obs.offset_indices = offset_indices


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

    # Patch missing parameters if possible, otherwise break:
    if atm.tpars is None and ret.itemp is not None:
        if len(ret.map_pars['temp']) < atm.temp_model.npars:
            log.error('Not all temp parameters are defined (tpars)')
        atm.tpars = np.zeros(atm.temp_model.npars)
        atm.tpars[ret.map_pars['temp']] = ret.params[ret.itemp]
    if atm.molpars is None and ret.imol is not None:
        if len(ret.map_pars['mol']) < len(atm.mol_pnames):
            log.error('Not all abundance parameters are defined (molpars)')
        atm.molpars = np.zeros(len(atm.mol_pnames))
        atm.molpars[ret.map_pars['mol']] = ret.params[ret.imol]
    if pyrat.rayleigh.pars is None and ret.iray is not None:
        if len(ret.map_pars['ray']) < pyrat.rayleigh.npars:
            log.error('Not all Rayleigh parameters are defined (rpars)')
        pyrat.rayleigh.pars = np.zeros(pyrat.rayleigh.npars)
        pyrat.rayleigh.pars[ret.map_pars['ray']] = ret.params[ret.iray]
    if pyrat.cloud.pars is None and ret.icloud is not None:
        if len(ret.map_pars['cloud']) < pyrat.cloud.npars:
            log.error('Not all Cloud parameters are defined (cpars)')
        pyrat.cloud.pars = np.zeros(pyrat.cloud.npars)
        pyrat.cloud.pars[ret.map_pars['cloud']] = ret.params[ret.icloud]

