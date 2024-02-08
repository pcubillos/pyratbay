# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Retrieval',
]

import numpy as np

from .. import constants as pc
from .. import tools as pt


class ColorTheme():
    """
    Descriptor object that remembers the default input theme
    (in obj._default_theme) until it is redefined.
    Used in pyrat.plot_spectrum().

    To understand this sorcery see:
    https://docs.python.org/3/howto/descriptor.html
    """
    def __set_name__(self, obj, name):
        self.private_name = '_' + name

    def __get__(self, obj, objtype=None):
        value = getattr(obj, self.private_name)
        return value

    def __set__(self, obj, value):
        priv_name = self.private_name
        setattr(obj, priv_name, value)
        obj._default_theme = value


class Retrieval():
    theme = ColorTheme()

    def __init__(self, inputs, atm, phy, obs, opacity, log):
        self.nparams = 0     # Number of free parameters
        self.posterior = None
        self.bestp = None
        self.spec_best = None
        self.spec_low1 = None
        self.spec_low2 = None
        self.spec_high1 = None
        self.spec_high2 = None

        self.mcmcfile = inputs.mcmcfile
        self.retflag = inputs.retflag
        self.qcap = inputs.qcap
        if atm.chemistry == 'tea':
            self.qcap = None
        # Lower and upper temperature retrieval boundaries
        self.tlow = inputs.tlow
        self.thigh = inputs.thigh

        self.sampler = inputs.sampler
        theme = 'blue' if inputs.theme is None else inputs.theme
        self.theme = pt.resolve_theme(theme)
        # If defaulted to None, keep _default_theme as None until
        # the user redefines ret.theme to something else:
        if inputs.theme is None:
            self._default_theme = None
        # Retrieval configuration
        self.statistics = inputs.statistics
        self.nsamples = inputs.nsamples
        self.burnin = inputs.burnin
        self.thinning = inputs.thinning
        self.nchains = inputs.nchains
        self.grbreak = inputs.grbreak
        self.grnmin = inputs.grnmin
        self.resume = inputs.resume
        self.nlive = inputs.nlive

        # Overrides retflag. At some point this will be the only way.
        if inputs.retrieval_params is not None:
            pars = [
                par for par in inputs.retrieval_params.splitlines()
                if par != ''
            ]
            nparams = len(pars)
            pnames = []
            params = np.zeros(nparams)
            pmin = np.tile(-np.inf, nparams)
            pmax = np.tile(np.inf, nparams)
            pstep = np.zeros(nparams)
            prior = np.zeros(nparams)
            priorlow = np.zeros(nparams)
            priorup = np.zeros(nparams)
            for i,par in enumerate(pars):
                fields = par.split()
                nfields = len(fields)
                if nfields not in [2, 5, 7, 8]:
                    log.error(
                        "Invalid number of fields for retrieval_params entry\n"
                        f"'{par}'"
                    )
                pnames.append(fields[0])
                params[i] = fields[1]
                if nfields == 2:
                    continue
                pmin[i] = fields[2]
                pmax[i] = fields[3]
                pstep[i] = fields[4]
                if nfields == 5:
                    continue
                prior[i] = fields[5]
                priorlow[i] = fields[6]
                if nfields == 7:
                    priorup[i] = fields[6]
                else:
                    priorup[i] = fields[7]
            self.pnames = pnames
            self.params = params
            self.pstep = pstep
            self.pmin = pmin
            self.pmax = pmax
            self.prior = prior
            self.priorlow = priorlow
            self.priorup = priorup
        else:
            self.params = inputs.params
            self.pstep = inputs.pstep
            self.pmin = inputs.pmin
            self.pmax = inputs.pmax
            self.prior = inputs.prior
            self.priorlow = inputs.priorlow
            self.priorup = inputs.priorup
            self.pnames = collect_pnames_from_retflag(
                inputs, self, atm, opacity, log,
            )

        unique_pnames, counts = np.unique(self.pnames, return_counts=True)
        if np.any(counts>1):
            log.error(f'Repeated parameter names: {unique_pnames[counts>1]}')

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

        temp_pnames = []
        if atm.tmodelname in pc.tmodels:
            temp_pnames = atm.temp_model.pnames

        opacity.is_patchy = (
            opacity.fpatchy is not None or
            'f_patchy' in self.pnames
        )

        opacity_pnames = []
        for names in opacity.pnames:
            opacity_pnames += names

        offset_pnames = obs.offset_inst
        error_pnames = obs.uncert_scaling

        # Indices to map free parameters of each model:
        self.map_pars = map_pars = {
            'temp': [],
            'mol': [],
            'opacity': [[] for model in opacity.models],
            'offset': [],
            'error': [],
        }
        # Model parameter names
        self.nparams = len(self.pnames)
        self.texnames = [None for _ in self.pnames]

        solo_params = [
            'log_p_ref',
            'R_planet',
            'M_planet',
            'f_patchy',
            'T_eff',
            'f_dilution',
        ]
        all_available_params = (
            solo_params +
            temp_pnames +
            atm.mol_pnames +
            opacity_pnames +
            offset_pnames +
            error_pnames
        )

        # Indices for each model parameters in self.params array:
        self.irad = None
        self.ipress = None
        self.ipatchy = None
        self.imass = None
        self.itstar = None
        self.idilut = None
        self.itemp = None
        self.imol = None
        self.iopacity = [[] for model in opacity.models]
        self.ioffset = None
        self.ierror = None

        itemp = []
        imol = []
        ioffset = []
        ierror = []
        for i,pname in enumerate(self.pnames):
            if pname == 'log_p_ref':
                self.ipress = np.array([i])
                self.texnames[i] = r'$\log p_{{\rm ref}}$'
            elif pname == 'R_planet':
                self.irad = np.array([i])
                self.texnames[i] = fr'$R_{{\rm p}}$ ({utex[atm.runits]})'
            elif pname == 'M_planet':
                self.imass = np.array([i])
                self.texnames[i] = fr'$M_{{\rm p}}$ ({utex[atm.mass_units]})'
            elif pname =='f_patchy':
                self.ipatchy = np.array([i])
                self.texnames[i] = r'$\phi_{\rm patchy}$'
            elif pname == 'T_eff':
                self.itstar = np.array([i])
                self.texnames[i] = r'$T_{\rm eff}$ (K)'
            elif pname == 'f_dilution':
                self.idilut = np.array([i])
                self.texnames[i] = r'$F_{\rm dilut}$'

            elif pname in temp_pnames:
                itemp.append(i)
                idx = atm.temp_model.pnames.index(pname)
                map_pars['temp'].append(idx)
                self.texnames[i] = atm.temp_model.texnames[idx]
            elif pname in atm.mol_pnames:
                imol.append(i)
                for imodel,model in enumerate(atm.vmr_models):
                    if pname in model.pnames:
                        idx = model.pnames.index(pname)
                        map_pars['mol'].append((imodel,idx))
                        self.texnames[i] = atm.mol_texnames[imodel]
                        break
            elif pname in opacity_pnames:
                for j,model in enumerate(opacity.models):
                    if pname in opacity.pnames[j]:
                        self.iopacity[j].append(i)
                        idx = model.pnames.index(pname)
                        map_pars['opacity'][j].append(idx)
                        self.texnames[i] = model.texnames[idx]
                        if not np.isfinite(model.pars[idx]):
                            model.pars[idx] = self.params[i]
                        break
            elif pname in offset_pnames:
                ioffset.append(i)
                idx = offset_pnames.index(pname)
                map_pars['offset'].append(idx)
                self.texnames[i] = obs.depth.offset_texnames[idx]
            elif pname in error_pnames:
                ierror.append(i)
                idx = error_pnames.index(pname)
                map_pars['error'].append(idx)
                self.texnames[i] = obs.depth.err_texnames[idx]
            else:
                log.error(
                    f"Invalid retrieval parameter '{pname}'. Possible "
                    f"values are:\n{all_available_params}"
                )

        if len(itemp) > 0:
            self.itemp = itemp
        if len(imol) > 0:
            self.imol = imol
        if len(ioffset) > 0:
            self.ioffset = ioffset
        if len(ierror) > 0:
            self.ierror = ierror

        # Patch missing parameters if possible:
        patch_temp = (
            atm.tpars is None and
            self.itemp is not None and
            len(map_pars['temp']) == atm.temp_model.npars
        )
        if patch_temp:
            atm.tpars = np.zeros(atm.temp_model.npars)
            atm.tpars[map_pars['temp']] = self.params[self.itemp]


        patch_abundance = (
            atm.vmr_pars is None and
            self.imol is not None and
            len(map_pars['mol']) == atm.mol_npars
        )
        if patch_abundance:
            atm.vmr_pars = [
                [None for _ in range(vmr_model.npars)]
                for vmr_model in atm.vmr_models
            ]
            for j,imol in enumerate(self.imol):
                imodel,idx = map_pars['mol'][j]
                atm.vmr_pars[imodel][idx] = params[imol]

        if self.ipatchy is not None and opacity.fpatchy is None:
            opacity.fpatchy = self.params[self.ipatchy]

        if atm.tpars is None and atm.temp_model is not None:
            log.error('Not all temperature parameters were defined (tpars)')

        if atm.vmr_pars is None and atm.mol_npars > 0:
            log.error('Not all vmr parameter values were defined (vmr_vars)')
        bad_models = ''
        for j,model in enumerate(opacity.models):
            if hasattr(model, 'pars') and not np.all(np.isfinite(model.pars)):
                bad_models = f"{opacity.models_type[j]} model '{model.name}', "
        if len(bad_models) > 0:
            log.error(f'Undefined parameter values for {bad_models[:-2]}')


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('Retrieval information:')
        if self.params is None:
            fw.write('No retrieval parameters set.')
            return fw.text

        pmin = [None for _ in self.params] if self.pmin is None else self.pmin
        pmax = [None for _ in self.params] if self.pmax is None else self.pmax
        psteps = [None for _ in self.params] if self.pstep is None else self.pstep

        fw.write(
            '  Parameter name        value        pmin        pmax       pstep',
        )
        fw.write(
            '  {:15}  {:>10}  {:>10}  {:>10}  {:>10}',
            '(pnames)', '(params)', '(pmin)', '(pmax)', '(pstep)',
        )
        n_obs = len(self.pnames)

        for i in range(n_obs):
            fw.write(
                '  {:15s}  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}',
                self.pnames[i], self.params[i], pmin[i], pmax[i], psteps[i],
            )

        if self.prior is not None:
            fw.write('\nParameter name     Prior')
            for i in range(n_obs):
                pname = self.pnames[i]
                val = self.params[i]
                if psteps[i] == 0.0:
                    fw.write('  {:15s}  Fixed at  {:10.3e}', pname, val)
                elif psteps[i] < 0.0:
                    j = -int(psteps[i]) - 1
                    fw.write('  {:15s}  Shared with  {}', pname, self.pnames[j])
                elif self.priorlow[i]==0 and self.priorup[i]==0:
                    fw.write(
                        '  {:15s}  Uniform between     [{:10.3e}, {:10.3e}]',
                        pname, pmin[i], pmax[i],
                    )
                else:
                    fw.write(
                        '  {:15s}  Gaussian  {:10.3e} -{:.3e}  {:+.3e}',
                        pname, self.prior[i], self.priorlow[i], self.priorup[i],
                    )

        fw.write('\nRetrieval algorithm (sampler): {}', self.sampler)
        if self.sampler is None:
            return fw.text
        fw.write('Number of retrieval samples (nsamples): {:,}', self.nsamples)
        # if self.sampler == 'snooker':
        fw.write('Number of parallel chains (nchains):   {}', self.nchains)
        fw.write('Number of burned-in samples (burnin):  {:,}', self.burnin)
        fw.write('Thinning factor (thinning): {}', self.thinning)
        if self.grbreak > 0.0:
            fw.write(
                'Gelman-Rubin convergence criterion to stop (grbreak): '
                f'{self.grbreak}',
            )
            if self.grnmin > 1:
                fw.write(
                    'Minimum number of samples before GR stop (grnmin): '
                    f'{int(self.grnmin)}'
                )
            else:
                fw.write(
                    'Minimum fraction of samples before GR stop (grnmin): '
                    f'{self.grnmin}'
                )

        fw.write(
            f'\nUpper boundary for sum of metal abundances (qcap): {self.qcap}',
        )
        fw.write(f'Temperature upper boundary (tlow, K):  {self.tlow:6.1f}')
        fw.write(f'Temperature lower boundary (thigh, K): {self.thigh:6.1f}')

        fw.write('\nRetrieval posterior file (mcmcfile): {}', self.mcmcfile)
        if self.posterior is not None:
            nsamples, nparams = np.shape(self.posterior)
            fw.write(
                '\nParameter name     Best-fit   Posterior distribution '
                f'of shape [{nsamples},{nparams}]\n'
                '                   (bestp)    (posterior)',
            )
            post = iter(self.posterior.T)
            for pname, bestp, pstep in zip(self.pnames, self.bestp, psteps):
                if pstep > 0:
                    fw.write(
                        '  {:15} {:10.3e}  {}', pname, bestp, next(post),
                        fmt={'float':'{: .3e}'.format}, edge=2,
                    )
                else:
                    fw.write('  {:15} {:10.3e}', pname, bestp)
            fw.write(
                '\nBest-fit spectrum (spec_best):\n    {}', self.spec_best,
                fmt={'float':'{: .3e}'.format},
            )
        return fw.text


def collect_pnames_from_retflag(
        inputs, ret, atm, opacity, log,
    ):
    """
    Check and setup retrieval parameters via retflag argument.
    To be deprecated, use retrieval_parameters instead.
    """
    retflag = ret.retflag
    pnames = []

    if retflag is None:
        return pnames

    rayleigh_pnames = []
    for model in opacity.models:
        if inputs.rayleigh is None:
            continue
        if model.name in inputs.rayleigh and hasattr(model, 'pnames'):
            rayleigh_pnames += model.pnames

    cloud_pnames = []
    for model in opacity.models:
        if inputs.clouds is None:
            continue
        if model.name in inputs.clouds and hasattr(model, 'pnames'):
            cloud_pnames += model.pnames


    if 'temp' in retflag and atm.tmodelname is None:
        log.error('Requested temp in retflag, but there is no tmodel')
    if 'mol' in retflag:
        if atm.vmr_vars == []:
            log.error("Requested mol in retflag, but there is no 'vmr_vars'")
        if atm.bulk is None:
            log.error('Requested mol in retflag, but there are no bulk species')
    if 'ray' in retflag and inputs.rayleigh is None:
        log.error('Requested ray in retflag, but there are no rayleigh models')

    if 'cloud' in retflag and inputs.clouds is None:
        log.error('Requested cloud in retflag, but there are no cloud models')

    if 'temp' in retflag:
        pnames += atm.temp_model.pnames
    if 'rad' in retflag:
        pnames += ['R_planet']
    if 'press' in retflag:
        pnames += ['log_p_ref']
    if 'mol' in retflag:
        pnames += atm.mol_pnames
    if 'ray' in retflag:
        pnames += rayleigh_pnames
    if 'cloud' in retflag:
        pnames += cloud_pnames
    if 'patchy' in retflag:
        pnames += ['f_patchy']
    if 'mass' in retflag:
        pnames += ['M_planet']
    if 'tstar' in retflag:
        pnames += ['T_eff']
    if 'offset' in retflag:
        pnames += list(ret.offset_inst)
    if 'error' in retflag:
        pnames += list(ret.uncert_scaling)

    # Retrieval variables:
    if ret.params is not None and len(ret.params) != len(pnames):
        nparams = len(ret.params)
        log.error(
            f'The number of input fitting parameters (params, {nparams}) does '
            f'not match the number of required parameters ({len(pnames)})'
        )
    return pnames
