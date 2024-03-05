# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Hydrogen_Ion',
]

import numpy as np
from .. import constants as pc


class Hydrogen_Ion():
    """Hydrogen ion opacity model based on John (1988)"""
    def __init__(self, wn):
        """
        Parameters
        ----------
        wn: 1D float array
            Wavenumber array where to sample the opacities (cm-1).
        """
        self.name = 'H- bound-free/free-free'
        self.wn = wn
        self.species = ['H', 'e-']
        self.nwave = len(self.wn)
        self._alpha = pc.h * pc.c / pc.k
        # Photo-detachment wn threshold in cm-1 (wl0 = 1.6419 um):
        self._wn0_bf = 6090.5

        self.sigma_bf = self._sigma_bound_free()
        self._setup_free_free()


    def _sigma_bound_free(self):
        """
        Compute the bound-free cross section for H- in cm2 units.
        Follows Equation (4) of John 1988, AA, 193, 189.

        Examples
        --------
        >>> import pyratbay.opacity as op
        >>> import pyratbay.spectrum as ps
        >>> import pyratbay.constants as pc

        >>> wl = ps.constant_resolution_spectrum(0.1, 14.0, 15000.0)
        >>> h_ion = op.Hydrogen_Ion(1e4/wl)
        >>> # H- cross sections in cm2 / H-_molec and cm5 / elec / H_molec
        >>> temperature = 2000.0
        >>> sigma_bf = h_ion.cross_section_bound_free(temperature)
        >>> sigma_ff = h_ion.cross_section_free_free(temperature)

        >>> # As in MacDonald & Lewis (2022), Fig (6):
        >>> plt.figure('MacDonald & Lewis (2022)')
        >>> plt.clf()
        >>> ax = plt.subplot(111)
        >>> plt.plot(wl, sigma_bf, color='blue', label='H- bound free')
        >>> plt.plot(wl, sigma_ff, color='orange', label='H- free free')
        >>> plt.xscale('log')
        >>> plt.yscale('log')
        >>> ax.tick_params(which='both', right=True, top=True, direction='in')
        >>> ax.legend(loc='upper right')
        >>> plt.xlim(0.4, 14.0)
        >>> plt.ylim(1.0e-40, 1.0e-35)
        >>> plt.xlabel('Wavelength (um)')
        >>> plt.ylabel(r'Cross section (cm$^5$ / H / e-)')

        >>> # Lothringer (2018), Fig (6):
        >>> temp = np.array([4150.0, 3350.0, 2850.0])
        >>> press = np.array([1e-4, 1e-2, 1.0])

        >>> vmr_Hm = np.array([3.0e-12, 4.0e-11, 8e-11])
        >>> vmr_H = np.array([1.0, 1.0, 8e-2])
        >>> vmr_e = np.array([1.0e-4, 4.0e-6, 5.0e-7])
        >>> n_Hm = vmr_Hm * press*pc.bar / pc.k / temp
        >>> n_H = vmr_H * press*pc.bar / pc.k / temp
        >>> n_e = vmr_e * press*pc.bar / pc.k / temp
        >>> sigma_ff = h_ion.cross_section_free_free(temp)
        >>> dens = np.array([n_H,n_e]).T
        >>> extinction = h_ion.calc_extinction_coefficient(temp, dens)

        >>> cols = 'blue xkcd:green orange'.split()
        >>> labels = ['0.1 mbar','10 mbar', '1000 mbar']
        >>> plt.figure('Lothringer et al. (2018)')
        >>> plt.clf()
        >>> ax = plt.subplot(111)
        >>> for i in range(3):
        >>>     plt.plot(wl, extinction[i], c=cols[i], lw=1.5, label=labels[i])
        >>> plt.xscale('log')
        >>> plt.yscale('log')
        >>> ax.tick_params(which='both', right=True, top=True, direction='in')
        >>> plt.xlim(0.1, 5.0)
        >>> plt.ylim(1e-15, 1.0e-5)
        >>> ax.legend(loc='upper left')
        >>> plt.xlabel('Wavelength (um)')
        >>> plt.ylabel('Extinction coefficient (cm-1)')
        """
        # Bound-free photo-detachment cross section constants:
        c_bf = [152.519, 49.534, -118.858, 92.536, -34.194, 4.982]

        wn_mask = self.wn > self._wn0_bf
        reduced_wl = 1e-2*np.sqrt(self.wn[wn_mask] - self._wn0_bf)
        # Equation (5):
        f_lambda = np.zeros(np.sum(wn_mask))
        for n in range(6):
            f_lambda += c_bf[n] * reduced_wl**n

        # bound-free cross section in cm2 per H- molecule, Equation (4):
        sigma = np.zeros(len(self.wn))
        sigma[wn_mask] = 1.0e-6 * (reduced_wl/self.wn[wn_mask])**3.0 * f_lambda

        return sigma


    def _setup_free_free(self):
        """
        Equation (6) of John 1988, AA, 193, 189.

        Pre-compute this on initialization.
        """
        wl = 1e4/self.wn
        a_short = [ 518.1021,   473.2636, -482.2089, 115.5291]
        b_short = [-734.8666,  1443.4137, -737.1616, 169.6374]
        c_short = [1021.1775, -1977.3395, 1096.8827, -245.649]
        d_short = [-479.0721,   922.3575, -521.1341,  114.243]
        e_short = [  93.1373,  -178.9275,  101.7963, -21.9972]
        f_short = [  -6.4285,    12.3600,   -7.0571,   1.5097]

        a_long = [ 2483.346,  -3449.889,   2200.040,   -696.271,    88.283]
        b_long = [  285.827,  -1158.382,   2427.719,  -1841.400,   444.517]
        c_long = [-2054.291,   8746.523, -13651.105,   8624.970, -1863.864]
        d_long = [ 2827.776, -11485.632,  16755.524, -10051.530,  2095.288]
        e_long = [-1341.537,   5303.609,  -7510.494,   4400.067,  -901.788]
        f_long = [  208.952,   -812.939,   1132.738,   -655.020,   132.985]

        wl_crit = 0.3645
        sw_mask = (0.182 < wl) & (wl < wl_crit)
        lw_mask = wl >= wl_crit
        wl_short = wl[sw_mask]
        wl_long = wl[lw_mask]

        nwave_sw = np.sum(sw_mask)
        ff_sw_factors = np.zeros((4,nwave_sw))
        for i in range(4):
            ff_sw_factors[i] = (
                a_short[i] * wl_short**2.0 +
                b_short[i] +
                c_short[i] / wl_short +
                d_short[i] / wl_short**2.0 +
                e_short[i] / wl_short**3.0 +
                f_short[i] / wl_short**4.0
            )
        nwave_lw = np.sum(lw_mask)
        ff_lw_factors = np.zeros((5,nwave_lw))
        for i in range(5):
            ff_lw_factors[i] = (
                a_long[i]*wl_long**2.0 +
                b_long[i] +
                c_long[i]/wl_long +
                d_long[i]/wl_long**2.0 +
                e_long[i]/wl_long**3.0 +
                f_long[i]/wl_long**4.0
            )
        self._ff_sw_factors = 1.0e-29 * np.expand_dims(ff_sw_factors.T,axis=2)
        self._ff_lw_factors = 1.0e-29 * np.expand_dims(ff_lw_factors.T,axis=2)
        self._sw_mask = sw_mask
        self._lw_mask = lw_mask


    def cross_section_free_free(self, temperature):
        """
        Compute the free-free cross section for H- in cm5/H_mol/electron.

        Equation (6) of John 1988, AA, 193, 189.
        """
        scalar_temp = np.isscalar(temperature)
        if scalar_temp:
            temperature = np.array([temperature])
        # Do not go beyond valid ranges of polynomial domain:
        temperature = np.clip(temperature, 1000.0, 10080.0)

        sigma_ff = np.zeros((self.nwave,len(temperature)))
        beta = np.array([
            np.sqrt(5040.0/temperature)**(i+2)
            for i in range(6)
        ])

        sigma_ff[self._sw_mask] = np.sum(beta[0:4]*self._ff_sw_factors, axis=1)
        sigma_ff[self._lw_mask] = np.sum(beta[1:6]*self._ff_lw_factors, axis=1)
        # Cross section in cm5 / H_molecule / electron
        # (electron number density = e_pressure / kT)
        sigma_ff *= pc.k * temperature
        cross_section_ff = sigma_ff.T
        if scalar_temp:
            return cross_section_ff[0]
        return cross_section_ff


    def cross_section_bound_free(self, temperature):
        """
        Compute the bound-free cross section for H- in cm5/H_mol/electron.

        Equation (3) of John 1988, AA, 193, 189.
        """
        scalar_temp = np.isscalar(temperature)
        if not scalar_temp:
            temperature = np.expand_dims(temperature, axis=1)

        # Cross section in cm5 / H_molecule / electron
        cross_section_bf = (
            0.75 * temperature**-1.5 * pc.k
            * np.exp(self._wn0_bf * self._alpha/temperature)
            * (1.0-np.exp(-self.wn*self._alpha/temperature))
            * self.sigma_bf
        )

        return cross_section_bf


    def calc_extinction_coefficient(self, temperature, density, layer=None):
        """
        Calculate extinction coefficient spectra (cm-1) for given
        temperature and number-density profiles.

        Parameters
        ----------
        temperature: Float or 1D float array
            Temperature profile in Kelvin (of size ntemp)
        density: 1D/2D float array
            Number density profiles (molecules per cm3) of H and electron
            over the given temperature array.
            If temperature is 1D array, must be of shape [ntemp,2]
            If temperature is scalar, density must be of size 2.
        layer: Integer
            If not None, compute extinction coefficient only at selected layer.
            In this case, the output array dimension is reduced by 1.

        Returns
        -------
        extinction_coefficient: 1D/2D float array
            H- extinction coefficient spectrum (cm-1)
            If temperature is scalar, output is a 1D array of length nwave.
            Otherwise, outout is a 2D array of shape [ntemp,nwave]
        """
        # Dimentional checks:
        is_scalar = np.isscalar(temperature)
        if is_scalar:
            if np.size(density) != 2:
                raise ValueError(
                    'Incompatible dimensions, if temperature is scalar '
                    'density must have 2 elements (H and e- densities)'
                )
        elif np.shape(density) != (len(temperature), 2):
            ntemp = len(temperature)
            raise ValueError(
                'Incompatible dimensions, density must be a 2D array '
                f'of shape [{ntemp}, 2], i.e., [ntemp, 2]'
            )
        if not is_scalar and layer is not None:
            temperature = temperature[layer]
            density = density[layer]
            is_scalar = True

        cross_section = (
            self.cross_section_bound_free(temperature) +
            self.cross_section_free_free(temperature)
        )

        if is_scalar:
            return cross_section * np.prod(density)

        # Total H- extinction coefficient (cm-1)
        extinction_coefficient = (
            cross_section *
            np.prod(density, axis=1, keepdims=True)
        )
        return extinction_coefficient

