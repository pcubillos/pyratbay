# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'pressure',
    'temperature',
    'chemistry',
    'uniform',
    'hydro_g',
    'hydro_m',
    'hill_radius',
    'stoich',
    'mean_weight',
    'ideal_gas_density',
    'equilibrium_temp',
    'transit_path',
    'temperature_posterior',
]

import numpy as np
import scipy.integrate as si
import scipy.constants as sc
import scipy.interpolate as sip
import mc3.utils as mu
import chemcat as cat

from .. import tools as pt
from .. import constants as pc
from .. import io as io
from . import tmodels


def pressure(ptop, pbottom, nlayers, units="bar", log=None, verb=0):
    """
    Compute a log-scale pressure profile.

    Parameters
    ----------
    ptop: String or Float
        Pressure at the top of the atmosphere. If string, may contain units.
    pbottom: String or Float
       Pressure at the bottom of the atmosphere. If string, may contain units.
    nlayers: Integer
       Number of pressure layers.
    units: String
       Pressure input units (if not defined in ptop, pbottom).
       Available units are: barye, mbar, pascal, bar (default), and atm.
    log: Log object
       Screen-output log handler.
    verb: Integer
       Verbosity level (when log is None). Print out when verb > 0.

    Returns
    -------
    press: 1D float ndarray
       The pressure profile (in bars).

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants as pc

    >>> nlayers = 9
    >>> # These are all equivalent:
    >>> p1 = pa.pressure(ptop=1e-6, pbottom=1e2, nlayers=nlayers)
    >>> p2 = pa.pressure(1e-6, 1e2, nlayers, units='bar')
    >>> p3 = pa.pressure('1e-6 bar', '1e2 bar', nlayers)
    >>> p4 = pa.pressure(1e-6*pc.bar, 1e2*pc.bar, nlayers, units='barye')
    >>> print(p1)
    [1.e-06 1.e-05 1.e-04 1.e-03 1.e-02 1.e-01 1.e+00 1.e+01 1.e+02]
    """
    if log is None:
        log = mu.Log(verb=verb)
    # Unpack pressure input variables:
    ptop = pt.get_param(ptop, units, gt=0.0)
    pbottom = pt.get_param(pbottom, units, gt=0.0)

    ptop_bar = ptop/pc.bar
    pbottom_bar = pbottom/pc.bar

    if ptop >= pbottom:
        log.error(
            f'Bottom-layer pressure ({pbottom_bar:.2e} bar) must be '
            f'higher than the top-layer pressure ({ptop_bar:.2e} bar).'
        )

    # Create pressure array in bars:
    press = np.logspace(np.log10(ptop_bar), np.log10(pbottom_bar), nlayers)
    log.head(
        f'Creating {nlayers}-layer atmospheric model between '
        f'{pbottom_bar:.1e} and {ptop_bar:.1e} {units}.'
    )
    return press


def temperature(tmodel, pressure=None, nlayers=None, log=None, params=None):
    """
    Temperature profile wrapper.

    Parameters
    ----------
    tmodel: String
        Name of the temperature model.
    pressure: 1D float ndarray
        Atmospheric pressure profile in bars.
    nlayers: Integer
        Number of pressure layers.
    log: Log object
        Screen-output log handler.
    params: 1D float ndarray
        Temperature model parameters. If None, return a tuple with the
        temperature model, its arguments, and the number or required parameters.

    Returns
    -------
    If params is not None:
        temperature: 1D float ndarray
            The evaluated atmospheric temperature profile.
    If params is None:
        temp_model: Callable
            The atmospheric temperature model.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa

    >>> nlayers = 11
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, units="bar")
    >>> # Isothermal profile:
    >>> temp_iso = pa.temperature("isothermal", pressure, params=1500.0)
    >>> print(temp_iso)
    [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]

    >>> # Guillot (2010) temperature profile:
    >>> params = np.array([-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0])
    >>> temp = pa.temperature('guillot', pressure, params=params)
    >>> print(temp)
    [1046.89057381 1046.89090056 1046.89433798 1046.93040895 1047.30779086
     1051.21739055 1088.76131307 1312.57904127 1640.18896334 1659.78818839
     1665.09706555]
    """
    if log is None:
        log = mu.Log()

    if pressure is not None and nlayers is None:
        nlayers = len(pressure)
    if tmodel in pc.tmodels:
        temp_model = tmodels.get_model(tmodel, pressure=pressure)
    else:
        log.error(
            f"Invalid temperature model '{tmodel}'.  Select from: {pc.tmodels}"
        )

    if params is None:
        return temp_model
    else:
        if np.size(params) != temp_model.npars:
            log.error(
                f"Wrong number of parameters ({np.size(params)}) for the "
                f"{temp_model.name} temperature model ({temp_model.npars})."
            )
        temperature = temp_model(params)
        log.head(f'\nComputed {tmodel} temperature model.')
        return temperature


def uniform(abundances, nlayers):
    """
    Generate an atmospheric file with uniform abundances.
    Save it into atmfile.

    Parameters
    ----------
    abundances: 1D float ndarray
        The species mole mixing ratio.
    nlayers: Integer
       Number of pressure layers.

    Returns
    -------
    vmr: 2D Float ndarray
        Abundance profiles of shape [nlayers,nspecies]

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> nlayers = 11
    >>> species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    >>> abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    >>> vmr = pa.uniform(abundances, nlayers)
    >>> print(vmr)
    [[8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]]
    """
    # Expand abundances to 2D array:
    if nlayers < 1:
        raise ValueError('The number of layers has to be larger than zero')
    vmr = np.tile(abundances, (nlayers,1))
    return vmr


def chemistry(
        chem_model, pressure, temperature, species,
        metallicity=0.0, e_abundances={}, e_scale={}, e_ratio={},
        q_uniform=None,
        solar_file=None, log=None, verb=1,
        atmfile=None, punits='bar',
    ):
    """
    Compute atmospheric abundaces for given pressure and
    temperature profiles with either uniform abundances or TEA.

    Parameters
    ----------
    chem_model: String
        Name of chemistry model, select from: 'uniform' or 'tea'
    pressure: 1D float ndarray
        Atmospheric pressure profile (bars).
    temperature: 1D float ndarray
        Atmospheric temperature profile (Kelvin).
    species: 1D string list
        Output atmospheric composition.
    metallicity: Float
        Metallicity enhancement factor in dex units relative to solar.
    e_abundances: Dictionary
        Custom elemental abundances.
        The dict contains the name of the element and their custom
        abundance in dex units relative to H=12.0.
        These values override metallicity.
    e_scale: Dictionary
        Scaling abundance factor for specified atoms by the respective
        values (in dex units, in addition to metallicity scaling).
        E.g. (3x solar): e_scale = {'C': np.log10(3.0)}
    e_ratio: Dictionary
        Custom elemental abundances scaled relative to another element.
        The dict contains the pair of elements joined by an underscore
        and their ratio in dex units, e.g., for a C/O ratio of 0.8 set
        e_ratio = {'C_O': 0.8}.
        These values modify the abundances after metallicity and
        e_scale have been applied.
    solar_file: String
        Input solar elemental abundances file (default Asplund et al. 2021).
    log: Log object
        Screen-output log handler.
    verb: Integer
        Verbosity level.
    atmfile: String
        If not None, output file where to save the atmospheric model.
    punits: String
        Output pressure units.

    Returns
    -------
    model: Callable
        The atmospheric chemistry-network model.

    Example
    -------
    >>> import pyratbay.atmosphere as pa
    >>> import matplotlib.pyplot as plt
    >>> import pyratbay.plots as pp

    >>> nlayers = 100
    >>> T0 = 1500.0
    >>> pressure = pa.pressure('1e-8 bar', '1e3 bar', nlayers)
    >>> temperature = pa.temperature('isothermal', pressure, nlayers, params=T0)
    >>> species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He H+ e-'.split()
    >>> # Equilibrium abundances model:
    >>> chem_model = 'tea'
    >>> chem_network = pa.chemistry(chem_model, pressure, temperature, species)

    >>> q_uniform = np.array([
    >>>     5e-4, 3e-5, 2e-4, 1e-8,  1e-6, 1e-14, 1e-13, 5e-10,
    >>>     1e-4, 0.85, 5e-3, 0.14,  3e-23, 1e-23])
    >>> chem_model = 'uniform'
    >>> chem_network_unif = pa.chemistry(
    >>>     chem_model, pressure, temperature, species, q_uniform=q_uniform)

    >>> # Plot the results:
    >>> ax1 = pp.abundance(
    >>>     chem_network.vmr, pressure, chem_network.species,
    >>>     colors='default', xlim=[1e-30, 3.0])

    >>> ax2 = pp.abundance(
    >>>     chem_network_unif.vmr, pressure, chem_network_unif.species,
    >>>     colors='default', xlim=[1e-30, 3.0], fignum=506)
    """
    if solar_file is None:
        solar_file = 'asplund_2021'
    if log is None:
        log = mu.Log(verb=verb)

    # Safety first
    nlayers = len(pressure)
    if len(temperature) != nlayers:
        log.error(
            f"pressure ({nlayers}) and temperature array "
            f"lengths ({len(temperature)}) don't match"
        )


    log.head("\nCompute chemical abundances.")
    chem_network = cat.Network(
        pressure, temperature, species,
        metallicity=metallicity,
        e_abundances=e_abundances,
        e_scale=e_scale,
        e_ratio=e_ratio,
        e_source=solar_file,
    )

    if chem_model == 'uniform':
        if len(species) != len(q_uniform):
            log.error(
                f"Species ({len(species)}) and q_uniform "
                f"array lengths ({len(q_uniform)}) don't match"
            )
        abundances = [
            q_uniform[list(species).index(spec)]
            for spec in chem_network.species
        ]
        chem_network.vmr = uniform(abundances, nlayers)


    elif chem_model == 'tea':
        chem_network.thermochemical_equilibrium()

    if atmfile is not None:
        header = "# TEA atmospheric file\n\n"
        io.write_atm(
            atmfile, pressure, temperature,
            chem_network.species, chem_network.vmr,
            punits=punits, header=header,
        )

    return chem_network


def hydro_g(pressure, temperature, mu, g, p0=None, r0=None):
    """
    Calculate radii using the hydrostatic-equilibrium equation considering
    a constant gravity.

    Parameters
    ----------
    pressure: 1D float ndarray
        Atmospheric pressure for each layer (in bar).
    temperature: 1D float ndarray
        Atmospheric temperature for each layer (in K).
    mu: 1D float ndarray
        Mean molecular mass for each layer (in g mol-1).
    g: Float
        Atmospheric gravity (in cm s-2).
    p0: Float
        Reference pressure level (in bar) where radius(p0) = r0.
    r0: Float
        Reference radius level (in cm) corresponding to p0.

    Returns
    -------
    radius: 1D float ndarray
        Radius for each layer (in cm).

    Notes
    -----
    If the reference values (p0 and r0) are not given, set radius = 0.0
    at the bottom of the atmosphere.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants as pc
    >>> nlayers = 11
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    >>> temperature = pa.tmodels.Isothermal(pressure)(1500.0)
    >>> mu = np.tile(2.3, nlayers)
    >>> g = pc.G * pc.mjup / pc.rjup**2
    >>> r0 = 1.0 * pc.rjup
    >>> p0 = 1.0  # bar
    >>> # Radius profile in Jupiter radii:
    >>> radius = pa.hydro_g(pressure, temperature, mu, g, p0, r0) / pc.rjup
    >>> print(radius)
    [1.0563673  1.04932138 1.04227547 1.03522956 1.02818365 1.02113774
     1.01409182 1.00704591 1.         0.99295409 0.98590818]
    """
    # Apply the HE equation:
    radius = si.cumulative_trapezoid(
        -pc.k*sc.N_A*temperature / (mu*g),
        np.log(pressure)
    )
    radius = np.concatenate(([0.0], radius))

    # Set absolute radii values if p0 and r0 are provided:
    if p0 is not None and r0 is not None:
        # Find current radius at p0:
        radinterp = sip.interp1d(pressure, radius, kind='slinear')
        r0_interp = radinterp(p0)
        # Set: radius(p0) = r0
        radius += r0 - r0_interp
    # Set radius = 0 at the bottom of the atmosphere:
    else:
        radius -= radius[-1]

    return radius


def hydro_m(pressure, temperature, mu, mass, p0, r0):
    """
    Calculate radii using the hydrostatic-equilibrium equation considering
    a variable gravity: g(r) = G*mass/r**2

    Parameters
    ----------
    pressure: 1D float ndarray
        Atmospheric pressure for each layer (in bar).
    temperature: 1D float ndarray
        Atmospheric temperature for each layer (in K).
    mu: 1D float ndarray
        Mean molecular mass for each layer (in g mol-1).
    mass: Float
        Object's mass (in g).
    p0: Float
        Reference pressure level (in bar) where radius(p0) = r0.
    r0: Float
        Reference radius level (in cm) corresponding to p0.

    Returns
    -------
    radius: 1D float ndarray
        Radius for each layer (in cm).

    Notes
    -----
    It is possible that this hydrostatic solution diverges when an
    atmosphere is too puffy.  In such cases, some returned radii
    will have np.inf values (at the top layers).

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants  as pc

    >>> nlayers = 11
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    >>> temperature = pa.tmodels.Isothermal(pressure)(1500.0)
    >>> mu = np.tile(2.3, nlayers)
    >>> mplanet = 1.0 * pc.mjup
    >>> r0 = 1.0 * pc.rjup
    >>> p0 = 1.0  # bar
    >>> # Radius profile in Jupiter radii:
    >>> radius = pa.hydro_m(pressure, temperature, mu, mplanet, p0, r0)/pc.rjup
    >>> print(radius)
    [1.05973436 1.05188019 1.04414158 1.036516   1.029001   1.02159419
     1.01429324 1.00709591 1.         0.99300339 0.986104  ]
    """
    # Apply the HE equation:
    I = si.cumulative_trapezoid(
        pc.k*sc.N_A*temperature/(pc.G*mu*mass),
        np.log(pressure)
    )
    I = np.concatenate(([0.0], I))

    # Find current radius at p0:
    radinterp = sip.interp1d(pressure, I, kind='slinear')
    I0 = radinterp(p0)
    # Set: radius(p0) = r0
    radius = 1.0/(I - I0 + 1/r0)

    # Divergent profile (radius should be monotonically decreasing):
    for j in range(len(radius)-1):
        if radius[j] <= radius[j+1]:
            radius[0:j+1] = np.inf

    return radius


def hill_radius(smaxis, mplanet, mstar):
    """
    Compute the Hill radius.  If any argument is None, return inf.

    Parameters
    ----------
    smaxis: Float
        Orbital semi-major axis (in cm).
    mplanet: Float
        Planetary mass (in g).
    mstar: Float
        Stellar mass (in g).

    Returns
    -------
    rhill: Float
        Hill radius of planet.
    """
    if smaxis is None or mplanet is None or mstar is None:
        return np.inf
    rhill = smaxis * (mplanet/(3.0*mstar))**(1.0/3.0)
    return rhill


def stoich(species):
    """
    Extract the elemental composition from a list of species.

    Parameters
    ----------
    species: 1D string list
        List of species.

    Returns
    -------
    elements: 1D string list
        List of elements contained in species list.
    stoich: 2D integer ndarray
        Stoichiometric elemental values for each species (number of elements).

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> species = ['H2', 'He', 'H2O', 'CO', 'CO2', 'CH4']
    >>> elements, stoichs = pa.stoich(species)
    >>> print(f'{elements}\n{stoichs}')
    ['C', 'H', 'He', 'O']
    [[0 2 0 0]
     [0 0 1 0]
     [0 2 0 1]
     [1 0 0 1]
     [1 0 0 2]
     [1 4 0 0]]
    """
    # Elemental composition and quantity for each species:
    comp, n = [], []

    for spec in species:
        comp.append([])
        n.append([])
        for char in spec:
            # Special case for electrons:
            if spec == 'e-':
                comp[-1].append(spec)
                n[-1].append(1)
            # New element:
            elif char.isupper():
                comp[-1].append(char)
                n[-1].append(1)
            # Same element:
            elif char.islower():
                comp[-1][-1] += char
            # Quantity:
            elif char.isdigit():
                n[-1][-1] = int(char)

    # Flatten nested list, and get unique set of elements:
    elements = sorted(set([
        element
        for spec in comp
        for element in spec
    ]))
    stoich = np.zeros((len(species), len(elements)), int)

    # Count how many elements in each species:
    for i in range(len(species)):
        for j in range(len(comp[i])):
            idx = elements.index(comp[i][j])
            stoich[i,idx] = n[i][j]

    return elements, stoich


def mean_weight(abundances, species=None, molfile=None, mass=None):
    """
    Calculate the mean molecular weight (a.k.a. mean molecular mass)
    for the given abundances composition.

    Parameters
    ----------
    abundances: 2D float iterable
        Species volume mixing fraction, of shape [nlayers,nmol].
    species: 1D string iterable
        Species names.
    molfile: String
        A molecules file with the species info.  If None, use
        pyratbay's default molecules.dat file.
    mass: 1D float ndarray
        The mass for each one of the species (g mol-1).
        If not None, this variable takes precedence over molfile.

    Returns
    -------
    mu: 1D float ndarray
        Mean molecular weight at each layer for the input abundances.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> species = 'H2 He H2O CO CO2 CH4'.split()
    >>> abundances = [[0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]]
    >>> mu = pa.mean_weight(abundances, species)
    >>> print(mu)
    [2.31939114]
    """
    if mass is None and species is None:
        raise ValueError('Either species or mass arguments must be specified')
    if mass is None:
        if molfile is None:
            molfile = pc.ROOT + 'pyratbay/data/molecules.dat'
        names, mass, radius = io.read_molecs(molfile)

        if np.any(~np.isin(species, names)):
            with pt.log_error():
                missing = species[~np.isin(species, names)]
                raise ValueError(f"Missing species masses for {missing}")
        mass = np.array([mass[names==spec][0] for spec in species])

    return np.sum(np.atleast_2d(abundances)*mass, axis=1)


def ideal_gas_density(abundances, pressure, temperature):
    """
    Use the ideal gas law to calculate number density in molecules cm-3.

    Parameters
    ----------
    abundances: 2D float array
        Species volume mixing fraction.
        Can have a 2D shape of [nlayers,nmol] or a 1D shape of [nlayers]
    pressure: 1D array
        Atmospheric pressure profile (in bars).
    temperature: 1D array
        Atmospheric temperature profile (in kelvin).

    Returns
    -------
    density: 2D float array
        Atmospheric density in molecules cm-3.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> nlayers = 11
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    >>> temperature = np.tile(1500.0, nlayers)
    >>> species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    >>> abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    >>> vmr = pa.uniform(abundances, nlayers)
    >>> dens = pa.ideal_gas_density(vmr, pressure, temperature)
    >>> print(dens[0])
    [4.10241993e+10 7.24297303e+09 4.82864869e+06 4.82864869e+06
     4.82864869e+02 4.82864869e+06]
    """
    if np.shape(abundances) == np.shape(pressure):
        return abundances * (pressure*pc.bar) / (temperature*pc.k)
    return abundances*np.expand_dims(pressure/temperature,axis=1) * pc.bar/pc.k


def equilibrium_temp(
        tstar, rstar, smaxis,
        A=0.0, f=1.0,
        tstar_unc=0.0, rstar_unc=0.0, smaxis_unc=0.0,
    ):
    r"""
    Calculate equilibrium temperature and uncertainty.

    Parameters
    ----------
    tstar: Scalar
        Effective temperature of host star (in kelvin degrees).
    rstar: Scalar
        Radius of host star (in cm).
    smaxis: Scalar
        Orbital semi-major axis (in cm).
    A: Scalar
        Planetary bond albedo.
    f: Scalar
        Planetary energy redistribution factor:
        f=0.5  no redistribution (total dayside reemission)
        f=1.0  good redistribution (4pi reemission)
    tstar_unc: Scalar
        Effective temperature uncertainty (in kelvin degrees).
    rstar_unc: Scalar
        Stellar radius uncertainty (in cm).
    smaxis_unc: Scalar
        Semi-major axis uncertainty (in cm).

    Returns
    -------
    teq: 1D ndarray
        Planet equilibrium temperature in kelvin degrees.
    teq_unc: 1D ndarray
        Equilibrium temperature uncertainty.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants as pc
    >>> import numpy as np
    >>> # HD 209458b (Stassun et al. 2017):
    >>> tstar  = 6091.0
    >>> rstar  = 1.19 * pc.rsun
    >>> smaxis = 0.04747 * pc.au
    >>> tstar_unc  = 10.0
    >>> rstar_unc  = 0.02 * pc.rsun
    >>> smaxis_unc = 0.00046 * pc.au
    >>> A = 0.3
    >>> f = np.array([1.0, 0.5])
    >>> teq, teq_unc = pa.equilibrium_temp(tstar, rstar, smaxis, A, f,
    >>>     tstar_unc, rstar_unc, smaxis_unc)
    >>> print(f'HD 209458b T_eq =\n    '
    >>>    f'{teq[0]:.1f} +/- {teq_unc[0]:.1f} K (day--night redistributed)\n'
    >>>    f'    {teq[1]:.1f} +/- {teq_unc[1]:.1f} K (instant re-emission)')
    HD 209458b T_eq =
        1345.1 +/- 13.2 K (day--night redistribution)
        1599.6 +/- 15.7 K (instant re-emission)
    """
    teq = ((1.0-A)/f)**0.25 * (0.5*rstar/smaxis)**0.5 * tstar

    teq_unc = teq * np.sqrt(
        (    tstar_unc/tstar)**2.0
      + (0.5*smaxis_unc/smaxis)**2.0
      + (0.5*rstar_unc/rstar)**2.0
    )

    return teq, teq_unc


def transit_path(radius, nskip=0):
    """
    Calculate the distances between layers for a set of rays grazing
    an atmosphere defined by concentric spheres at radii radius.
    Assume that the grazing rays have an impact parameter = radius.
    See for example, Fig 1 of Molliere et al. (2019), AA, 627, 67.

    Parameters
    ----------
    radius: 1D float ndarray
        Atmospheric radius profile array (from top to bottom).
    nskip: Integer
        Number of layers to skip from the top of the atmosphere.

    Returns
    -------
    path: List of 1D float ndarrays
        List where each element is a 1D array of the paths at an
        impact parameter defined.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> nlayers = 5
    >>> radius = np.linspace(5.0, 1.0, nlayers)
    >>> path = pa.transit_path(radius)
    >>> # First path grazes the outer layer (thus, empty path)
    >>> # Second path is sqrt(5.0**2 - 4.0**2), and so on
    >>> for p in path:
    >>>     print(p)
    []
    [3.]
    [1.35424869 2.64575131]
    [1.11847408 1.22803364 2.23606798]
    [1.02599614 1.04455622 1.09637632 1.73205081]
    >>> # Now, ignore the top layer:
    >>> path = pa.transit_path(radius, nskip=1)
    >>> for p in path:
    >>>     print(p)
    []
    []
    [2.64575131]
    [1.22803364 2.23606798]
    [1.04455622 1.09637632 1.73205081]
    """
    rad = radius[nskip:]
    nlayers = len(rad)

    # Empty-filling layers that don't contribute:
    path = []
    for r in range(nskip):
        path.append(np.empty(0, np.double))

    # Compute the path for each impact parameter:
    r = 0
    while r < nlayers:
        raypath = np.empty(r, np.double)
        for i in range(r):
            raypath[i] = (
                np.sqrt(rad[i  ]**2 - rad[r]**2) -
                np.sqrt(rad[i+1]**2 - rad[r]**2)
            )
        path.append(raypath)
        r += 1

    return path


def temperature_posterior(posterior, tmodel):
    """
    Compute the median and inter-quantiles regions (68% and 95%)
    of a temperature profile posterior.

    Parameters
    ----------
    posterior: 2D float ndarray
        A posterior distribution for the parameters of tmodel
        of shape [nsamples, npars].
    tmodel: Callable
        Temperature-profile model.

    Returns
    -------
    median: 1D float ndarray [nlayers]
        The matplotlib Axes of the figure.
    low1: 1D float ndarray [nlayers]
        Lower temperature boundary of the 68%-interquantile range.
    high1: 1D float ndarray [nlayers]
        Upper temperature boundary of the 68%-interquantile range.
    low2: 1D float ndarray [nlayers]
        Lower temperature boundary of the 95%-interquantile range.
    high2: 1D float ndarray [nlayers]
        Upper temperature boundary of the 95%-interquantile range.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.plots as pp
    >>> import numpy as np

    >>> # Non-inverted temperature profile:
    >>> pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers=100)
    >>> tmodel = pa.tmodels.Guillot(pressure)

    >>> # Simulate posterior where log_kappa' and log_gamma1 vary:
    >>> nsamples = 5000
    >>> posterior = np.tile([-4.0, -1.0, 0.0, 0.0, 1000.0, 0.0], (nsamples,1))
    >>> posterior[:,0] = np.random.normal(-4.0, 0.5, nsamples)
    >>> posterior[:,1] = np.random.normal(-1.0, 0.5, nsamples)

    >>> # Compute posterior and take a look at it:
    >>> tpost = pa.temperature_posterior(posterior, tmodel)
    >>> ax = pp.temperature(pressure, profiles=tpost[0], bounds=tpost[1:])
    """
    u, uind, uinv = np.unique(
        posterior[:,0], return_index=True, return_inverse=True,
    )
    nlayers = len(tmodel.pressure)
    nsamples = len(u)

    # Evaluate posterior PT profiles:
    profiles = np.zeros((nsamples, nlayers), np.double)
    for i in range(nsamples):
        profiles[i] = tmodel(posterior[uind[i]])

    # Get median and inter-quantile ranges for 68% and 95% percentiles:
    median = np.zeros(nlayers, np.double)
    low1 = np.zeros(nlayers, np.double)
    low2 = np.zeros(nlayers, np.double)
    high1 = np.zeros(nlayers, np.double)
    high2 = np.zeros(nlayers, np.double)
    for i in range(nlayers):
        tpost = profiles[uinv,i]
        median[i] = np.percentile(tpost, 50.000)
        low2[i] = np.percentile(tpost,  2.275)
        low1[i] = np.percentile(tpost, 15.865)
        high1[i] = np.percentile(tpost, 84.135)
        high2[i] = np.percentile(tpost, 97.725)
    return median, low1, high1, low2, high2

