# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Hydrogen_Ion_Opacity',
]

import numpy as np
from .. import constants as pc


class Hydrogen_Ion_Opacity():
    """Hydrogen ion opacity model based on John (1988)"""
    def __init__(self, wl, species=None):
        """
        Parameters
        ----------
        wl: 1D float array
            Wavelength array where to compute the opacities (in microns).
        species: 1D string iterable
            List of atmospheric species names.  If provided, this will
            pre-load the indices of H, H-, and e-, which are needed
            later to compute extinction coefficients.
        """
        self.wl = wl
        self.ec = None
        self.cross_section_bf = self.bound_free_cross_section()
        self.free_free_setup()

        # Pre-load indices for a known atmosphere:
        self.has_opacity = False
        if species is not None:
            self.has_opacity = np.all(np.in1d(('H', 'H-', 'e-'), species))
            self._h_ion_idx = np.where(np.array(species) == 'H-')[0]
            self._h_idx = np.where(np.array(species) == 'H')[0]
            self._e_idx = np.where(np.array(species) == 'e-')[0]

    def bound_free_cross_section(self):
        """
        Compute the bound-free cross section for H- in cm2 units.
        Follows Equation (4) of John 1988, AA, 193, 189.

        Examples
        --------
        >>> import pyratbay.opacity.hydrogen_ion as h
        >>> import numpy as np

        >>> wl = np.logspace(-1.0, 2.0, 1000)
        >>> h_ion = h.Hydrogen_Ion_Opacity(wl)
        >>> # H- cross sections in cm2 / H-_molec and cm5 / elec / H_molec
	>>> temperature = 2000.0
        >>> sigma_bf = h_ion.cross_section_bf
        >>> sigma_ff = h_ion.free_free_cross_section(temperature)

        >>> # As in MacDonald & Lewis (2022), Fig (6):
        >>> import pyratbay.constants as pc
        >>> plt.figure('MacDonald & Lewis (2022)')
        >>> plt.clf()
        >>> ax = plt.subplot(111)
        >>> plt.plot(wl, sigma_bf, color='blue', label='H- bound free')
        >>> plt.plot(wl, sigma_ff[0], color='orange', label='H- free free')
        >>> plt.xscale('log')
        >>> plt.yscale('log')
        >>> ax.tick_params(which='both', right=True, top=True, direction='in')
        >>> ax.legend(loc='upper right')
        >>> plt.xlim(0.4, 14.0)
        >>> plt.ylim(1e-40, 1.0e-14)
        >>> plt.xlabel('Wavelength (um)')
        >>> plt.ylabel('Cross section (cm2/H- or cm5/H/e-)')

        >>> # Lothringer (2018), Fig (6):
	>>> temperature = np.array([4150.0, 3350.0, 2850.0])
        >>> press = np.array([1e-4, 1e-2, 1.0]) * pc.bar

        >>> vmr_Hm = np.array([3.0e-12, 4.0e-11, 8e-11])
        >>> vmr_H = np.array([1.0, 1.0, 8e-2])
        >>> vmr_e = np.array([1.0e-4, 4.0e-6, 5.0e-7])
        >>> n_Hm = vmr_Hm * press / pc.k / temperature
        >>> n_H = vmr_H * press / pc.k / temperature
        >>> n_e = vmr_e * press / pc.k / temperature
        >>> sigma_ff = h_ion.free_free_cross_section(temperature)

        >>> cols = 'blue xkcd:green orange'.split()
        >>> labels = ['0.1 mbar','10 mbar', '1000 mbar']
        >>> plt.figure('Lothringer et al. (2018)')
        >>> plt.clf()
        >>> ax = plt.subplot(111)
        >>> for i in range(3):
        >>>     sigma = sigma_bf*n_Hm[i] + (sigma_ff[i]*n_H[i]*n_e[i])
        >>>     plt.plot(wl, sigma, c=cols[i], lw=2, label=labels[i])
        >>> plt.xscale('log')
        >>> plt.yscale('log')
        >>> ax.tick_params(which='both', right=True, top=True, direction='in')
        >>> plt.xlim(0.1, 5.0)
        >>> plt.ylim(1e-15, 1.0e-5)
        >>> ax.legend(loc='upper left')
        >>> plt.xlabel('Wavelength (um)')
        >>> plt.ylabel('Extinction coefficient (cm-1)')
        """
        wl = self.wl
        # Photo-detachment threshold (microns):
        wl0_bf = 1.6419
        # Bound-free photo-detachment cross section constants:
        c_bf = [152.519, 49.534, -118.858, 92.536, -34.194, 4.982]

        wl_mask = wl < wl0_bf
        reduced_wl = np.sqrt(1.0/wl[wl_mask] - 1/wl0_bf)
        # Equation (5):
        f_lambda = np.zeros(np.sum(wl_mask))
        for n in range(6):
            f_lambda += c_bf[n] * reduced_wl**n

        # bound-free cross section in cm2 per H- molecule, Equation (4):
        sigma = np.zeros(len(wl))
        sigma[wl_mask] = 1.0e-18 * wl[wl_mask]**3.0 * reduced_wl**3.0 * f_lambda

        return sigma


    def free_free_setup(self):
        """
        Equation (6) of John 1988, AA, 193, 189.
        """
        wl = self.wl
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
        self.ff_sw_factors = 1.0e-29 * np.expand_dims(ff_sw_factors.T,axis=2)
        self.ff_lw_factors = 1.0e-29 * np.expand_dims(ff_lw_factors.T,axis=2)
        self.sw_mask = sw_mask
        self.lw_mask = lw_mask


    def free_free_cross_section(self, temperature):
        """
        Compute the free-free cross section for H- in cm5/H_mol/electron.

        Equation (6) of John 1988, AA, 193, 189.
        """
        if np.isscalar(temperature):
            temperature = np.array([temperature])
        # Do not go beyond valid ranges of polynomial domain:
        temperature = np.clip(temperature, 1000.0, 10080.0)

        sigma_ff = np.zeros((len(self.wl),len(temperature)))
        beta = np.sqrt(5040.0/temperature)
        beta_arr = np.array([beta**(i+2) for i in range(6)])

        sigma_ff[self.sw_mask] = np.sum(beta_arr[0:4]*self.ff_sw_factors,axis=1)
        sigma_ff[self.lw_mask] = np.sum(beta_arr[1:6]*self.ff_lw_factors,axis=1)
        # Cross section in cm5 / H_molecule / electron
        # (but how? physics explanation TBD)
        sigma_ff *= pc.k * temperature
        self.sigma_ff = sigma_ff.T
        return self.sigma_ff


    def absorption(self, temperature, number_density, species=None):
        """
        Evaluate the H- extinction coefficient spectrum (cm-1)
        over the given atmospheric model.

        If either H, H-, or e- are not present in the species
        the H- extinction coefficient will not be computed and set
        to None.

        Parameters
        ----------
        temperature: 1D float array
            A temperature profile (in K).
        number_density: 2D float array
            The atmospheric number densities array (molecule cm-3)
            with shape nlayers (same as temperature) by nspecies
            (same as species list).
        species: 1D string iterable
            The list of atmospheric species.  If this was preset on
            initialization, it can be skipped here.
        """
        # Get species indices if needed:
        if species is not None:
            self.has_opacity = np.all(np.in1d(('H', 'H-', 'e-'), species))
            self._h_ion_idx = np.where(np.array(species) == 'H-')[0]
            self._h_idx = np.where(np.array(species) == 'H')[0]
            self._e_idx = np.where(np.array(species) == 'e-')[0]

        if not self.has_opacity:
            self.ec = None
            return

        # Total H- extinction coefficient (cm-1)
        ec_bf = self.cross_section_bf * number_density[:,self._h_ion_idx]
        ec_ff = (
            number_density[:,self._h_idx] *
            number_density[:,self._e_idx] *
            self.free_free_cross_section(temperature)
        )
        self.ec = ec_bf + ec_ff

