# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'pressure',
    'temperature',
    'uniform',
    'abundance',
    'hydro_g',
    'hydro_m',
    'rhill',
    'stoich',
    'mean_weight',
    'ideal_gas_density',
    'equilibrium_temp',
    'transit_path',
    'temperature_posterior',
    'make_atomic',
    'make_preatm',
    ]

import os
import subprocess

import numpy as np
import scipy.integrate as si
import scipy.constants as sc
import scipy.interpolate as sip
import mc3.utils as mu

from .. import tools as pt
from .. import constants as pc
from .. import io as io
from .  import tmodels


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
       The pressure profile (in barye units).

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants  as pc
    >>> nlayers = 9
    >>> # These are all equivalent:
    >>> p1 = pa.pressure(ptop=1e-6,   pbottom=1e2, nlayers=nlayers)
    >>> p2 = pa.pressure(1e-6,        1e2,         nlayers, 'bar')
    >>> p3 = pa.pressure('1e-6 bar', '1e2 bar',    nlayers)
    >>> p4 = pa.pressure(1e-6*pc.bar, 1e2*pc.bar,  nlayers, 'barye')
    >>> print(p1/pc.bar)
    [1.e-06 1.e-05 1.e-04 1.e-03 1.e-02 1.e-01 1.e+00 1.e+01 1.e+02]
    """
    if log is None:
        log = mu.Log(verb=verb)
    # Unpack pressure input variables:
    ptop    = pt.get_param(ptop,    units, gt=0.0)
    pbottom = pt.get_param(pbottom, units, gt=0.0)

    ptop_txt = ptop/pt.u(units)
    pbot_txt = pbottom/pt.u(units)

    if ptop >= pbottom:
        log.error(f'Bottom-layer pressure ({pbot_txt:.2e} {units}) must be '
            f'higher than the top-layer pressure ({ptop_txt:.2e} {units}).')

    # Create pressure array in barye (CGS) units:
    press = np.logspace(np.log10(ptop), np.log10(pbottom), nlayers)
    log.head(f'Creating {nlayers}-layer atmospheric model between '
             f'{pbot_txt:.1e} and {ptop_txt:.1e} {units}.')
    return press


def temperature(tmodel, pressure=None, nlayers=None, log=None, params=None):
    """
    Temperature profile wrapper.

    Parameters
    ----------
    tmodel: String
        Name of the temperature model.
    pressure: 1D float ndarray
        Atmospheric pressure profile in barye units.
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
    >>> # Isothermal profile:
    >>> temp_iso = pa.temperature("isothermal", params=1500.0, nlayers=nlayers)
    >>> print(temp_iso)
    [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]
    >>> # Three-channel Eddington-approximation profile:
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, "bar")
    >>> params = np.array([-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0])
    >>> temp = pa.temperature('tcea', pressure, params=params)
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
        temp_model = tmodels.get_model(
            tmodel, pressure=pressure, nlayers=nlayers)
    else:
        log.error(f"Invalid temperature model '{tmodel}'.  "
                  f"Select from: {pc.tmodels}")

    if params is None:
        return temp_model
    else:
        if np.size(params) != temp_model.npars:
            log.error(
                f"Wrong number of parameters ({np.size(params)}) for the "
                f"{temp_model.name} temperature model ({temp_model.npars}).")
        temperature = temp_model(params)
        log.head(f'\nComputed {tmodel} temperature model.')
        return temperature


def uniform(pressure, temperature, species, abundances, punits="bar",
            log=None, atmfile=None):
    """
    Generate an atmospheric file with uniform abundances.
    Save it into atmfile.

    Parameters
    ----------
    pressure: 1D float ndarray
        Monotonously decreasing pressure profile (in punits).
    temperature: 1D float ndarray
        Temperature profile for pressure layers (in Kelvin).
    species: 1D string ndarray
        List of atmospheric species.
    abundances: 1D float ndarray
        The species mole mixing ratio.
    punits:  String
       Pressure units.
    log: Log object
        Screen-output log handler.
    atmfile: String
        If not None, filename to save atmospheric model.

    Returns
    -------
    qprofiles: 2D Float ndarray
        Abundance profiles of shape [nlayers,nspecies]

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> nlayers = 11
    >>> punits = 'bar'
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, punits)
    >>> tmodel = pa.tmodels.Isothermal(nlayers)
    >>> species    = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    >>> abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    >>> qprofiles = pa.uniform(pressure, tmodel(1500.0), species,
    >>>     abundances=abundances, punits=punits)
    >>> print(qprofiles)
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
    if log is None:
        log = mu.Log()

    nlayers = len(pressure)
    # Safety checks:
    if len(temperature) != nlayers:
        log.error(f"Pressure array length ({nlayers}) and temperature array "
            f"length ({len(temperature)}) don't match.")
    if len(species) != len(abundances):
        log.error(f"Species array length ({len(species)}) and abundances "
            f"array length ({len(abundances)}) don't match.")

    # Expand abundances to 2D array:
    qprofiles = np.tile(abundances, (nlayers,1))

    if atmfile is not None:
        header = ("# This is an atmospheric file with pressure, temperature,\n"
                  "# and uniform mole mixing ratio profiles.\n\n")
        io.write_atm(atmfile, pressure, temperature, species, qprofiles,
            punits=punits, header=header)

    return qprofiles


def abundance(pressure, temperature, species, elements=None,
        quniform=None, atmfile=None, punits='bar', xsolar=1.0,
        escale={}, solar_file=None, log=None, verb=1, ncpu=1):
    """
    Compute atmospheric abundaces for given pressure and
    temperature profiles with either uniform abundances or TEA.

    Parameters
    ----------
    pressure: 1D float ndarray
        Atmospheric pressure profile (barye).
    temperature: 1D float ndarray
        Atmospheric temperature profile (Kelvin).
    species: 1D string list
        Output atmospheric composition.
    elements: 1D strings list
        Input elemental composition (default to minimum list of elements
        required to form species).
    quniform: 1D float ndarray
        If not None, the output species abundances (isobaric).
    atmfile: String
        If not None, output file where to save the atmospheric model.
    punits: String
        Output pressure units.
    xsolar: Float
        Metallicity enhancement factor.
    escale: Dict
        Multiplication factor for specified atoms (dict's keys)
        by the respective values (on top of the xsolar scaling).
    solar_file: String
        Input solar elemental abundances file (Default Asplund et al. 2009).
    log: Log object
        Screen-output log handler.
    verb: Integer
        Verbosity level.
    ncpu: Integer
        Number of parallel CPUs to use in TEA calculation.

    Returns
    -------
    q: 2D float ndarray
       Atmospheric volume mixing fraction abundances of shape
       [nlayers, nspecies].

    Example
    -------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants  as pc
    >>> atmfile = "pbtea.atm"
    >>> nlayers = 100
    >>> press = np.logspace(-8, 3, nlayers) * pc.bar
    >>> temp  = 900+500/(1+np.exp(-(np.log10(press)+1.5)*1.5))
    >>> species = ['H2O', 'CH4', 'CO', 'CO2', 'NH3', 'C2H2', 'C2H4', 'HCN',
    >>>            'N2', 'H2', 'H', 'He']
    >>> # Automatically get 'elements' necessary from the list of species:
    >>> Q = pa.abundance(press, temp, species)
    """
    if solar_file is None:
        solar_file = pc.ROOT + "pyratbay/data/AsplundEtal2009.txt"
    if log is None:
        log = mu.Log(verb=verb)
    # Uniform-abundances profile:
    if quniform is not None:
        log.head("\nCompute uniform-abundances profile.")
        q = uniform(
            pressure, temperature, species, quniform, punits, log, atmfile)
        return q

    # TEA abundances:
    log.head("\nCompute TEA thermochemical-equilibrium abundances profile.")
    # Prep up files:
    atomic_file, patm = "PBatomicfile.tea", "PBpreatm.tea"
    make_atomic(xsolar, escale, atomic_file, solar_file)
    if elements is None:
       elements, dummy = stoich(species)
    specs = elements + list(np.setdiff1d(species, elements))
    make_preatm(
        pressure/pt.u(punits), temperature, atomic_file, elements, specs, patm)
    # Run TEA:
    pt.make_tea(abun_file=atomic_file, verb=verb, ncpu=ncpu)
    proc = subprocess.Popen([pc.ROOT+"pyratbay/TEA/tea/runatm.py", patm, "TEA"])
    proc.communicate()

    # Reformat the TEA output into the pyrat format:
    if atmfile is None:
        atmfile = 'TEA.tea'
    io.import_tea("TEA.tea", atmfile, species)
    q = io.read_atm(atmfile)[4]
    os.remove(atomic_file)
    os.remove(patm)
    os.remove("TEA.cfg")
    os.remove("TEA.tea")
    return q


def hydro_g(pressure, temperature, mu, g, p0=None, r0=None):
    """
    Calculate radii using the hydrostatic-equilibrium equation considering
    a constant gravity.

    Parameters
    ----------
    pressure: 1D float ndarray
        Atmospheric pressure for each layer (in barye).
    temperature: 1D float ndarray
        Atmospheric temperature for each layer (in K).
    mu: 1D float ndarray
        Mean molecular mass for each layer (in g mol-1).
    g: Float
        Atmospheric gravity (in cm s-2).
    p0: Float
        Reference pressure level (in barye) where radius(p0) = r0.
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
    >>> temperature = pa.tmodels.Isothermal(nlayers)(1500.0)
    >>> mu = np.tile(2.3, nlayers)
    >>> g = pc.G * pc.mjup / pc.rjup**2
    >>> r0 = 1.0 * pc.rjup
    >>> p0 = 1.0 * pc.bar
    >>> # Radius profile in Jupiter radii:
    >>> radius = pa.hydro_g(pressure, temperature, mu, g, p0, r0) / pc.rjup
    >>> print(radius)
    [1.0563673  1.04932138 1.04227547 1.03522956 1.02818365 1.02113774
     1.01409182 1.00704591 1.         0.99295409 0.98590818]
    """
    # Apply the HE equation:
    radius = si.cumtrapz(-pc.k*sc.N_A*temperature / (mu*g), np.log(pressure))
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
        Atmospheric pressure for each layer (in barye).
    temperature: 1D float ndarray
        Atmospheric temperature for each layer (in K).
    mu: 1D float ndarray
        Mean molecular mass for each layer (in g mol-1).
    mass: Float
        Object's mass (in g).
    p0: Float
        Reference pressure level (in barye) where radius(p0) = r0.
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
    >>> temperature = pa.tmodels.Isothermal(nlayers)(1500.0)
    >>> mu = np.tile(2.3, nlayers)
    >>> mplanet = 1.0 * pc.mjup
    >>> r0 = 1.0 * pc.rjup
    >>> p0 = 1.0 * pc.bar
    >>> # Radius profile in Jupiter radii:
    >>> radius = pa.hydro_m(pressure, temperature, mu, mplanet, p0, r0)/pc.rjup
    >>> print(radius)
    [1.05973436 1.05188019 1.04414158 1.036516   1.029001   1.02159419
     1.01429324 1.00709591 1.         0.99300339 0.986104  ]
    """
    # Apply the HE equation:
    I = si.cumtrapz((pc.k*sc.N_A*temperature)/(pc.G*mu*mass), np.log(pressure))
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


def rhill(smaxis, mplanet, mstar):
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
    rhill = smaxis * (mplanet/(3*mstar))**(1.0/3.0)
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
    >>> print('{}\n{}'.format(elements, stoichs))
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
    elements = sorted(set([element for spec    in comp
                                   for element in spec]))
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
    >>> species     = ['H2', 'He', 'H2O', 'CO', 'CO2', 'CH4']
    >>> abundances  = [[0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]]
    >>> mu = pa.mean_weight(abundances, species)
    >>> print(mu)
    [2.31928918]
    """
    if mass is None and species is None:
        raise ValueError('Either species or mass arguments must be specified')
    if mass is None:
        if molfile is None:
            molfile = pc.ROOT + 'pyratbay/data/molecules.dat'
        names, mass, diam = io.read_molecs(molfile)
        mass = np.array([mass[names==spec][0] for spec in species])

    return np.sum(np.atleast_2d(abundances)*mass, axis=1)


def ideal_gas_density(abundances, pressure, temperature):
    """
    Use the Ideal gas law to calculate number density in molecules cm-3.

    Parameters
    ----------
    abundances: 2D float ndarray
        Species volume mixing fraction, of shape [nlayers,nmol].
    pressure: 1D ndarray
        Atmospheric pressure profile (in barye units).
    temperature: 1D ndarray
        Atmospheric temperature profile (in kelvin).

    Returns
    -------
    density: 2D float ndarray
        Atmospheric density in molecules cm-3.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> atmfile = "uniform_test.atm"
    >>> nlayers = 11
    >>> pressure    = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    >>> temperature = np.tile(1500.0, nlayers)
    >>> species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    >>> abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    >>> qprofiles = pa.uniform(pressure, temperature, species, abundances)
    >>> dens = pa.ideal_gas_density(qprofiles, pressure, temperature)
    >>> print(dens[0])
    [4.10241993e+10 7.24297303e+09 4.82864869e+06 4.82864869e+06
     4.82864869e+02 4.82864869e+06]
    """
    return abundances * np.expand_dims(pressure/temperature, axis=1) / pc.k


def equilibrium_temp(tstar, rstar, smaxis, A=0.0, f=1.0,
    tstar_unc=0.0, rstar_unc=0.0, smaxis_unc=0.0):
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
            raypath[i] = (np.sqrt(rad[i  ]**2 - rad[r]**2) -
                          np.sqrt(rad[i+1]**2 - rad[r]**2) )
        path.append(raypath)
        r += 1

    return path


def temperature_posterior(posterior, tmodel, tpars, ifree, pressure):
    """
    Compute the median and inter-quantiles regions (68% and 95%)
    of a temperature profile posterior.

    Parameters
    ----------
    posterior: 2D float ndarray [nsamples, nfree]
        A posterior distribution for tmodel.
    tmodel: Callable
        Temperature-profile model.
    tpars: 1D float iterable [npars]
        Temperature-profile parameters (including fixed parameters).
    ifree: 1D bool iterable [npars]
        Mask of free (True) and fixed (False) parameters in tpars.
        The number of free parameters must match nfree in posterior.
    pressure: 1D float ndarray [nlayers]
        The atmospheric pressure profile in barye.

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
    >>> import pyratbay.constants as pc
    >>> import pyratbay.plots as pp
    >>> import numpy as np

    >>> # Non-inverted temperature profile:
    >>> pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers=100)
    >>> tmodel = pa.tmodels.TCEA(pressure)
    >>> tpars = np.array([-4.0, -1.0, 0.0, 0.0, 1000.0, 0.0])
    >>> # Simulate posterior where kappa' and gamma are variable:
    >>> nsamples = 5000
    >>> ifree = [True, True, False, False, False, False]
    >>> posterior = np.array([
    >>>     np.random.normal(tpars[0], 0.5, nsamples),
    >>>     np.random.normal(tpars[1], 0.1, nsamples)]).T
    >>> tpost = pa.temperature_posterior(
    >>>     posterior, tmodel, tpars, ifree, pressure)
    >>> ax = pp.temperature(pressure, profiles=tpost[0], bounds=tpost[1:])
    """
    nlayers = len(pressure)

    u, uind, uinv = np.unique(
        posterior[:,0], return_index=True, return_inverse=True)
    nsamples = len(u)

    # Evaluate posterior PT profiles:
    profiles = np.zeros((nsamples, nlayers), np.double)
    tpars = np.array(tpars)
    ifree = np.array(ifree)
    for i in range(nsamples):
        tpars[ifree] = posterior[uind[i]]
        profiles[i] = tmodel(tpars)

    # Get median and inter-quantile ranges for 68% and 95% percentiles:
    low1   = np.zeros(nlayers, np.double)
    low2   = np.zeros(nlayers, np.double)
    median = np.zeros(nlayers, np.double)
    high1  = np.zeros(nlayers, np.double)
    high2  = np.zeros(nlayers, np.double)
    for i in range(nlayers):
        tpost = profiles[uinv,i]
        low2[i]   = np.percentile(tpost,  2.275)
        low1[i]   = np.percentile(tpost, 15.865)
        median[i] = np.percentile(tpost, 50.000)
        high1[i]  = np.percentile(tpost, 84.135)
        high2[i]  = np.percentile(tpost, 97.725)
    return median, low1, high1, low2, high2


def make_atomic(xsolar=1.0, escale={}, atomic_file=None, solar_file=None):
    """
    Generate an (atomic) elemental-abundances file by scaling a
    solar abundance composition.

    Parameters
    ----------
    xsolar: Integer
        Multiplication factor for metal abundances (everything
        except H and He).
    escale: Dict
        Multiplication factor for specified atoms (dict's keys)
        by the respective values (on top of the xsolar scaling).
    atomic_file: String
        Output filename of modified atomic abundances.
    solar_file: String
        Base (input) solar abundances filename. Defaulted to
        Solar atomic composition by Asplund et al. (2009)

    Returns
    -------
    z: 1D integer ndarray
        Atomic number
    symbol: 1D string ndarray
        Atom symbols
    dex: 1D float ndarray
        Logarithmic abundance respective to hydrogen, where log(H) = 12.
    name: 1D string ndarray
        Atom names.
    mass: 1D float ndarray
        Atomic mass in amu.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> # Compute and store elemental composition with [M/H] = 0.1:
    >>> afile = 'sub_solar_elemental_abundance.txt'
    >>> z, symbol, dex, names, mass = pa.make_atomic(
    >>>     xsolar=0.1, atomic_file=afile)
    >>> # Compute elemental composition with [M/H] = 10 and [C/H] = 100:
    >>> escale = {'C': 10.0}
    >>> z, symbol, dex, names, mass = pa.make_atomic(xsolar=10, escale=escale)
    """
    if solar_file is None:
        solar_file = pc.ROOT + 'pyratbay/data/AsplundEtal2009.txt'
    # Read the Asplund et al. (2009) solar elementa abundances:
    index, symbol, dex, name, mass = io.read_atomic(solar_file)

    # Scale the metals aundances:
    imetals = np.where((symbol != "H") & (symbol != "He"))
    dex[imetals] += np.log10(xsolar)

    for atom, fscale in escale.items():
        dex[symbol == atom] += np.log10(fscale)

    # Save data to file:
    if atomic_file is not None:
        with open(atomic_file, "w") as f:
            f.write("# Elemental abundances:\n"
                    "# atomic number, symbol, dex abundance, name, mass (u).\n")
            nelements = len(symbol)
            for i in range(nelements):
                f.write("{:3d}  {:2s}  {:5.2f}  {:10s}  {:12.8f}\n".
                    format(index[i], symbol[i], dex[i], name[i], mass[i]))
    return index, symbol, dex, name, mass


def make_preatm(pressure, temp, afile, elements, species, patm):
    """
    Generate a pre-atm file for TEA containing the elemental abundances
    at each atmospheric layer.

    Parameters
    ----------
    pressure: String
        Pressure atmospheric profile (bar).
    temp: 1D float array
        Temperature atmospheric profile (in K).
    afile: String
        Name of the elemental abundances file.
    elements: List of strings
        List of input elemental species.
    species: List of strings
        List of output molecular species.
    patm: String
        Output pre-atmospheric filename.

    Uncredited developers
    ---------------------
    Jasmina Blecic
    """
    # Number of layers:
    nlayers = len(pressure)

    # Get the elemental abundace data:
    index, symbol, dex, name, mass = io.read_atomic(afile)
    # Take only the elements we need:
    iatoms = np.in1d(symbol, elements)

    # Lisf of fractional number of molecules relative to hydrogen:
    nfrac = (10.0**(dex-dex[np.where(symbol=="H")]))[iatoms].tolist()

    # pyrat--TEA name dictionary:
    pyrat_to_tea = {}
    for line in open(pc.ROOT+"pyratbay/data/TEA_gdata_defaults.txt", "r"):
        pyrat_name, tea_name = line.split()
        pyrat_to_tea[pyrat_name] = tea_name

    # Check species names:
    for i in range(len(species)):
        if species[i].find("_") < 0:
            species[i] = pyrat_to_tea[species[i]]

    # Write pre-atm file:
    with open(patm, 'w') as f:
        # Pre-atm header with basic instructions
        f.write("# TEA pre-atmosphere file.\n"
                "# pressure (bar), temperature (K), abundance (unitless).\n\n")

        # Write species names:
        f.write(f'#SPECIES\n{" ".join(species)}\n\n')

        # Write TEA data:
        f.write("#TEADATA\n#Pressure          Temp  " +
                "  ".join([f"{atom:>12s}" for atom in symbol[iatoms]]) + "\n")
        for i in range(nlayers):
            f.write(f"{pressure[i]:10.4e}     {temp[i]:>8.2f}  ")
            f.write("  ".join([f"{abun:12.6e}" for abun in nfrac]) + "\n")
