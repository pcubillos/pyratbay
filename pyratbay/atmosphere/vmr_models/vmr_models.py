# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

"""
Admitedly, these classes are an overkill for what they do, however,
this setup really simplifies the call in pb.atmosphere.vmr_scale()
"""

__all__ = [
    'MetalEquil',
    'ScaleEquil',
    'RatioEquil',
    'IsoVMR',
    'ScaleVMR',
    'SlantVMR',
]

import functools
from collections.abc import Iterable

import numpy as np

from ... import constants as pc


def check_params(func):
    """Decorator to check that the number of model parameters is correct."""
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        self = args[0]
        params = kwargs['params'] if 'params' in kwargs else args[1]
        if np.size(params) != self.npars:
            raise ValueError(
                f'Number of parameters ({np.size(params)}) does '
                f'not match the required number of parameters ({self.npars}) '
                f'of the {self.name} model')
        return func(*args, **kwargs)
    return new_func


class MetalEquil():
    """Metallicity model for a thermochemical equilibrium model"""
    def __init__(self):
        self.name = 'metal_equil'
        self.pnames = ['[M/H]']
        self.texnames = ['[M/H]']
        self.npars = len(self.pnames)
        self.type = 'equil'

    def __str__(self):
        return (
            f'VMR model name: {self.name}\n'
            f'Number of parameters: {self.npars}\n'
            f'Parameters: {self.pnames}\n'
        )


class ScaleEquil():
    """Elemental abundance model for a thermochemical equilibrium model"""
    def __init__(self, name):
        """
        Parameters
        ----------
        name: String
            The element's name, the format must be: [X/H]
            where 'X' is the element name (e.g.: [C/H], [Na/H], ...).
        """
        self.name = 'scale_equil'
        self.element = name[1:-3]
        self.pnames = [name]
        self.texnames = [name]
        self.npars = len(self.pnames)
        self.type = 'equil'

    def __str__(self):
        return (
            f'VMR model name: {self.name}\n'
            f'Number of parameters: {self.npars}\n'
            f'Parameters: {self.pnames}'
        )


class RatioEquil():
    """Elemental abundance model for a thermochemical equilibrium model"""
    def __init__(self, name):
        """
        Parameters
        ----------
        name: String
            The elements name ratio, the format must be: X/Y
            where 'X' and 'Y' are the element names (e.g.: C/O, N/O, ...).
        """
        self.name = 'ratio_equil'
        idx = name.index('/')
        self.elements = [name[0:idx], name[idx+1:]]
        self.element_ratio = name.replace('/', '_')
        self.pnames = [name]
        self.texnames = [name]
        self.npars = len(self.pnames)
        self.type = 'equil'

    def __str__(self):
        return (
            f'VMR model name: {self.name}\n'
            f'Number of parameters: {self.npars}\n'
            f'Parameters: {self.pnames}'
        )


class IsoVMR():
    """Isobaric VMR model"""
    def __init__(self, species, pressure):
        """
        Parameters
        ----------
        species: String
            The atmospheric species name for this VMR profile model.
        pressure: 1D float iterable
            Pressure array (bar) where to evaluate the temperature profile.
        """
        self.species = species
        self.name = f'log_{species}'
        self.pnames = [f'log_{species}']
        self.texnames = [fr'$\log\ X_{{\rm {species}}}$']
        self.npars = len(self.pnames)
        self.pressure = pressure
        self.vmr = np.tile(1.0e-20, len(pressure))
        self.type = 'free'

    @check_params
    def __call__(self, params):
        """
        Parameters
        ----------
        params: scalar or iterable
            Desired isobaric log10(VMR).

        Returns
        -------
        vmr: 1D float ndarray
            VMR profile at each pressure.

        Examples
        --------
        >>> import pyratbay.atmosphere as pa

        >>> nlayers = 21
        >>> pressure = pa.pressure('1e-7 bar', '100 bar', nlayers)
        >>> iso_vmr_H2O = pa.vmr_models.IsoVMR('H2O', pressure)
        >>> # Evaluate model at VMR = 10**-3.0:
        >>> vmr_H2O = iso_vmr_H2O(-3.0)
        >>> print(vmr_H2O)
        [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
         0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001]
        """
        if isinstance(params, Iterable):
            self.vmr[:] = 10.0**params[0]
        else:
            self.vmr[:] = 10.0**params
        return np.copy(self.vmr)

    def __str__(self):
        return (
            f'VMR model name: {self.name}\n'
            f'Number of parameters: {self.npars}\n'
            f'Parameters: {self.pnames}'
        )


class ScaleVMR():
    """Scaled VMR model"""
    def __init__(self, species, pressure, vmr0):
        """
        Parameters
        ----------
        species: String
            The atmospheric species name for this VMR profile model.
        pressure: 1D float iterable
            Pressure array (bar) where to evaluate the temperature profile.
        vmr0: 1D float array
        """
        self.species = species
        self.name = f'scale_{species}'
        self.pnames = [f'log_{species}']
        self.texnames = [fr'$\log\ X_{{\rm {species}}}$']
        self.npars = len(self.pnames)
        self.pressure = pressure
        self.vmr0 = np.copy(vmr0)
        self.vmr = np.copy(self.vmr0)
        self.type = 'free'

    def __str__(self):
        return (
            f'VMR model name: {self.name}\n'
            f'Number of parameters: {self.npars}\n'
            f'Parameters: {self.pnames}'
        )

    @check_params
    def __call__(self, params):
        """
        Parameters
        ----------
        params: scalar or iterable
            Desired isobaric log10(VMR).

        Returns
        -------
        vmr: 1D float ndarray
            VMR profile at each pressure.

        Examples
        --------
        >>> import pyratbay.atmosphere as pa
        >>> import pyratbay.constants as pc

        >>> nlayers = 21
        >>> pressure = pa.pressure('1e-7 bar', '100 bar', nlayers)
        >>> # An initial VMR profile:
        >>> vmr_0 = 10**(0.5*np.tanh(np.log10(pressure)+2) - 3.5)

        >>> scale_vmr_H2O = pa.vmr_models.ScaleVMR('H2O', pressure, vmr_0)
        >>> # Reduce VMR by 1dex at each layer:
        >>> vmr_H2O = scale_vmr_H2O(-1.0)

        >>> plt.figure(0)
        >>> plt.clf()
        >>> plt.loglog(vmr_0, pressure)
        >>> plt.loglog(vmr_H2O, pressure)
        >>> plt.xlim(1e-6, 1e-2)
        >>> plt.ylim(100, 1e-7)
        """
        if isinstance(params, Iterable):
            self.vmr[:] = self.vmr0 * 10.0**params[0]
        else:
            self.vmr[:] = self.vmr0 * 10.0**params
        return np.copy(self.vmr)


class SlantVMR():
    """Slanted VMR model"""
    def __init__(self, species, pressure):
        """
        Parameters
        ----------
        species: String
            The atmospheric species name for this VMR profile model.
        pressure: 1D float iterable
            Pressure array (bar) where to evaluate the temperature profile.
        """
        self.species = species
        self.name = f'slant_{species}'
        self.pnames = [
            f'slope_{species}',
            f'log_VMR0_{species}',
            f'log_p0_{species}',
            f'min_log_{species}',
            f'max_log_{species}',
        ]
        self.texnames = [
            fr'$m_{{\rm {species}}}$',
            fr'$\log\ X_{{\rm {species}}}^{{0}}$',
            fr'$\log\ p_{{\rm {species}}}^{{0}}$',
            fr'$\log\ X_{{\rm {species}}}^{{\rm min}}$',
            fr'$\log\ X_{{\rm {species}}}^{{\rm max}}$',
        ]
        self.npars = len(self.pnames)
        self.pressure = pressure
        self.log_press = np.log10(self.pressure)
        self.vmr = np.tile(1.0e-20, len(pressure))
        self.type = 'free'

    def __str__(self):
        return (
            f'VMR model name: {self.name}\n'
            f'Number of parameters: {self.npars}\n'
            f'Parameters: {self.pnames}'
        )

    @check_params
    def __call__(self, params):
        """
        Parameters
        ----------
        params: 1D float ndarray
            VMR model parameters:
            slope:  VMR slope given as d_log(VMR) / d_log(p)
            vmr_0:  A reference log(VMR) value at log(p0)
            log_p0:  Reference pressure log(p0) where log(VMR)
            vmr_min:  Clip VMR profile a this minimum log(VMR) value
            vmr_max:  Clip VMR profile a this maximum log(VMR) value

        Returns
        -------
       vmr: 1D float ndarray
            VMR profile at each pressure.

        Examples
        --------
        >>> import pyratbay.atmosphere as pa
        >>> import pyratbay.constants as pc

        >>> nlayers = 51
        >>> pressure = pa.pressure('1e-7 bar', '100 bar', nlayers)
        >>> slant_vmr_H2O = pa.vmr_models.SlantVMR('H2O', pressure)

        >>> # Reduce VMR by 1 dex at each layer:
        >>> # slope, vmr0, log_p0, min_vmr, max_vmr
        >>> params = [0.0, -3.3, 0.0, -10, -1.0]
        >>> vmr_H2O = slant_vmr_H2O(params)
        >>> params = [0.5, -3.3, 0.0, -10, -1.0]
        >>> vmr_H2O_sloped = slant_vmr_H2O(params)
        >>> params = [0.75, -3.3, -3.0, -10, -1.0]
        >>> vmr_H2O_maxed = slant_vmr_H2O(params)

        >>> plt.figure(0)
        >>> plt.clf()
        >>> plt.loglog(vmr_H2O, pressure/pc.bar)
        >>> plt.loglog(vmr_H2O_sloped, pressure/pc.bar)
        >>> plt.loglog(vmr_H2O_maxed, pressure/pc.bar)
        >>> plt.xlim(1e-10, 1.0)
        >>> plt.ylim(100, 1e-7)
        """
        slope, vmr_0, log_p0, vmr_min, vmr_max = params
        log_vmr = slope * (self.log_press-log_p0) + vmr_0
        self.vmr[:] = 10.0**np.clip(log_vmr, vmr_min, vmr_max)
        return np.copy(self.vmr)

