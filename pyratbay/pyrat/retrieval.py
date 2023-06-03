# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Retrieval',
]

import numpy as np

from .. import tools as pt
from .. import constants as pc


class Retrieval():
    def __init__(self, inputs, atm, phy, obs, rayleigh, cloud, log):
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
        # Lower and upper temperature retrieval boundaries
        self.tlow = inputs.tlow
        self.thigh = inputs.thigh

        self.sampler = inputs.sampler
        self.theme = None
        # MCMC options
        self.nsamples = inputs.nsamples
        self.burnin = inputs.burnin
        self.thinning = inputs.thinning
        self.nchains = inputs.nchains
        self.grbreak = inputs.grbreak
        self.grnmin = inputs.grnmin
        self.resume = inputs.resume

        # Overrides retflag. At some point this will be the only way.
        if inputs.retrieval_params is not None:
            #self.pnames = []   # Model parameter names (screen)
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
            setup_retrieval_parameters_retflag(
                inputs, self, atm, obs, phy, rayleigh, cloud, log,
            )
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

        if atm.tmodelname in pc.tmodels:
            temp_pnames = atm.temp_model.pnames
        else:
            temp_pnames = []

        offset_pnames = obs.offset_instruments
        if offset_pnames is None:
            offset_pnames = []

        # Indices to map free parameters of each model:
        self.map_pars = map_pars = {
            'temp': [],
            'mol': [],
            'ray': [],
            'cloud': [],
            'offset': [],
        }
        # Model parameter names
        self.nparams = len(self.pnames)
        self.texnames = [None for _ in self.pnames]
        # Indices to parse the array of fitting parameters:

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
            rayleigh.pnames +
            cloud.pnames +
            offset_pnames
        )

        # Indices for each model parameters in self.params array:
        self.irad = None
        self.ipress = None
        self.ipatchy = None
        self.imass = None
        self.itstar = None
        self.itemp = None
        self.imol = None
        self.icloud = None
        self.iray = None
        self.ioffset = None

        itemp = []
        imol = []
        iray = []
        icloud = []
        ioffset = []
        for i,pname in enumerate(self.pnames):
            if pname == 'log_p_ref':
                self.ipress = np.array([i])
                self.texnames[i] = r'$\log p_{{\rm ref}}$'
            elif pname == 'R_planet':
                self.irad = np.array([i])
                self.texnames[i] = fr'$R_{{\rm p}}$ ({utex[atm.runits]})'
            elif pname == 'M_planet':
                self.imass = np.array([i])
                self.texnames[i] = fr'$M_{{\rm p}}$ ({utex[phy.mpunits]})'
            elif pname =='f_patchy':
                self.ipatchy = np.array([i])
                self.texnames[i] = r'$\phi_{\rm patchy}$'
            elif pname == 'T_eff':
                self.itstar = np.array([i])
                self.texnames[i] = r'$T_{\rm eff}$ (K)'

            elif pname in temp_pnames:
                itemp.append(i)
                idx = atm.temp_model.pnames.index(pname)
                map_pars['temp'].append(idx)
                self.texnames[i] = atm.temp_model.texnames[idx]
            elif pname in atm.mol_pnames:
                imol.append(i)
                idx = atm.mol_pnames.index(pname)
                map_pars['mol'].append(idx)
                self.texnames[i] = atm.mol_texnames[idx]
            elif pname in rayleigh.pnames:
                iray.append(i)
                idx = rayleigh.pnames.index(pname)
                map_pars['ray'].append(idx)
                self.texnames[i] = rayleigh.texnames[idx]
            elif pname in cloud.pnames:
                icloud.append(i)
                idx = cloud.pnames.index(pname)
                map_pars['cloud'].append(idx)
                self.texnames[i] = cloud.texnames[idx]
            elif pname in offset_pnames:
                ioffset.append(i)
                idx = offset_pnames.index(pname)
                map_pars['offset'].append(idx)
                self.texnames[i] = pname.replace('offset_', r'$\Delta$')
            else:
                log.error(
                    f"Invalid retrieval parameter '{pname}'. Possible "
                    f"values are:\n{all_available_params}"
                )
        # TBD: Check no repeated pnames

        if len(itemp) > 0:
            self.itemp = itemp
        if len(imol) > 0:
            self.imol = imol
        if len(icloud) > 0:
            self.icloud = icloud
        if len(iray) > 0:
            self.iray = iray
        if len(ioffset) > 0:
            self.ioffset = ioffset

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
            atm.molpars is None and
            self.imol is not None and
            len(map_pars['mol']) == atm.mol_npars
        )
        if patch_abundance:
            atm.molpars = np.zeros(len(atm.mol_pnames))
            atm.molpars[map_pars['mol']] = self.params[self.imol]

        patch_rayleigh = (
            rayleigh.pars is None and
            self.iray is not None and
            len(map_pars['ray']) == rayleigh.npars
        )
        if patch_rayleigh:
            rayleigh.pars = np.zeros(rayleigh.npars)
            rayleigh.pars[map_pars['ray']] = self.params[self.iray]

        patch_cloud = (
            cloud.pars is None and
            self.icloud is not None and
            len(map_pars['cloud']) == cloud.npars
        )
        if patch_cloud:
            cloud.pars = np.zeros(cloud.npars)
            cloud.pars[map_pars['cloud']] = self.params[self.icloud]

        if atm.tpars is None and atm.temp_model is not None:
            log.error('Not all temperature parameters were defined (tpars)')
        if atm.molpars is None and atm.mol_npars > 0:
            log.error('Not all abundance parameters were defined (molpars)')
        if cloud.pars is None and cloud.npars > 0:
            log.error('Not all Cloud parameters were defined (cpars)')
        if rayleigh.pars is None and rayleigh.npars > 0:
            log.error('Not all Rayleigh parameters were defined (rpars)')


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


def setup_retrieval_parameters_retflag(inputs, ret, atm, obs, phy, rayleigh, cloud, log):
    """
    Check and setup retrieval parameters via retflag argument.
    To be deprecated, use retrieval_parameters instead.
    """
    retflag = ret.retflag

    if retflag is None:
        return

    if 'temp' in retflag and atm.tmodelname is None:
        log.error('Requested temp in retflag, but there is no tmodel')
    if 'mol' in retflag:
        if inputs.molvars == []:
            log.error("Requested mol in retflag, but there is no 'molvars'")
        # TBD: This will break for pure eq-chem runs
        if atm.bulk is None:
            log.error('Requested mol in retflag, but there are no bulk species')
    if 'ray' in retflag and rayleigh.models == []:
        log.error('Requested ray in retflag, but there are no rayleigh models')
    if 'cloud' in retflag and cloud.models == []:
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
        ntemp = atm.temp_model.npars
        ret.itemp = np.arange(ret.nparams, ret.nparams + ntemp)
        ret.pnames += atm.temp_model.pnames
        ret.texnames += atm.temp_model.texnames
        ret.nparams += ntemp
        ret.map_pars['temp'] = np.arange(ntemp)
    if 'rad' in retflag:
        ret.irad = np.arange(ret.nparams, ret.nparams + 1)
        ret.pnames += ['R_planet']
        ret.texnames += [fr'$R_{{\rm p}}$ ({utex[atm.runits]})']
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
        nray = rayleigh.npars
        ret.iray = np.arange(ret.nparams, ret.nparams + nray)
        ret.pnames += rayleigh.pnames
        ret.texnames += rayleigh.texnames
        ret.nparams += nray
        ret.map_pars['ray'] = np.arange(nray)
    if 'cloud' in retflag:
        ncloud = cloud.npars
        ret.icloud = np.arange(ret.nparams, ret.nparams + ncloud)
        ret.pnames += cloud.pnames
        ret.texnames += cloud.texnames
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
        ret.texnames += [fr'$M_{{\rm p}}$ ({utex[phy.mpunits]})']
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

        band_names = [band.name for band in obs.filters]
        offset_indices = []
        for inst in ret.offset_instruments:
            flags = [inst in name for name in band_names]
            offset_indices.append(flags)
        obs.offset_indices = offset_indices


    # Retrieval variables:
    if ret.params is not None and len(ret.params) != ret.nparams:
        nparams = len(ret.params)
        log.error(
            f'The number of input fitting parameters (params, {nparams}) does '
            f'not match the number of required parameters ({ret.nparams})'
        )

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
    if rayleigh.pars is None and ret.iray is not None:
        if len(ret.map_pars['ray']) < rayleigh.npars:
            log.error('Not all Rayleigh parameters are defined (rpars)')
        rayleigh.pars = np.zeros(rayleigh.npars)
        rayleigh.pars[ret.map_pars['ray']] = ret.params[ret.iray]
    if cloud.pars is None and ret.icloud is not None:
        if len(ret.map_pars['cloud']) < cloud.npars:
            log.error('Not all Cloud parameters are defined (cpars)')
        cloud.pars = np.zeros(cloud.npars)
        cloud.pars[ret.map_pars['cloud']] = ret.params[ret.icloud]

