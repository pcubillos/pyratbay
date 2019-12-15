# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = [
    'makeatomic',
    'readatomic', 'makepreatm', 'TEA2pyrat', 'readmol',
    'pressure',
    'temperature',
    'uniform', 'abundances',
    'hydro_g', 'hydro_m',
    'stoich', 'meanweight', 'IGLdensity',
    'equilibrium_temp',
    ]

import os
import subprocess

import numpy as np
import scipy.integrate as si
import scipy.constants as sc
import scipy.interpolate as sip
import mc3.utils as mu

from .. import tools     as pt
from .. import constants as pc
from .. import io as io
from .  import tmodels


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
  >>> atmfile = "atm_test.dat"
  >>> nlayers = 11
  >>> punits  = 'bar'
  >>> pressure    = pa.pressure(1e-8, 1e2, nlayers, punits)
  >>> tmodel = pa.tmodels.Isothermal(nlayers)
  >>> species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
  >>> abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
  >>> qprofiles = pa.abundances(atmfile, pressure, tmodel(1500.0), species,
  >>>                           quniform=abundances, punits=punits)
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
      log.error("Pressure array length ({:d}) and Temperature array length "
                "({:d}) don't match.".format(nlayers, len(temperature)))
  if len(species) != len(abundances):
      log.error("Species array length ({:d}) and Abundances array length "
                "({:d}) don't match.".format(len(species), len(abundances)))

  # Expand abundances to 2D array:
  qprofiles = np.tile(abundances, (nlayers,1))

  if atmfile is not None:
      header = ("# This is an atmospheric file with pressure, temperature,\n"
                "# and uniform mole mixing ratio profiles.\n\n")
      io.write_atm(atmfile, pressure, temperature, species, qprofiles,
                   punits, header)

  return qprofiles


def makeatomic(solar, afile, xsolar=1.0, swap=None):
  """
  Generate an (atomic) elemental-abundances file by scaling the
  solar abundances from Asplund et al. (2009).
  http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A

  Parameters
  ----------
  solar: String
     Input solar abundances filename.
  afile: String
     Output filename of modified atomic abundances.
  xsolar: Integer
     Multiplication factor for metal abundances (everything
     except H and He).
  swap: 2-element string tuple
      Swap the abundances of the given elements.
  """
  # Read the Asplund et al. (2009) solar elementa abundances:
  index, symbol, dex, name, mass = readatomic(solar)
  # Count the number of elements:
  nelements = len(symbol)

  # Scale the metals aundances:
  imetals = np.where((symbol != "H") & (symbol != "He"))
  dex[imetals] += np.log10(xsolar)

  # Swap abundances if requested:
  if swap is not None:
      dex0 = dex[np.where(symbol == swap[0])]
      dex[np.where(symbol == swap[0])] = dex[np.where(symbol == swap[1])]
      dex[np.where(symbol == swap[1])] = dex0

  # Save data to file
  with open(afile, "w") as f:
      f.write("# Elemental abundances:\n"
              "# atomic number, symbol, dex abundance, name, molar mass (u).\n")
      for i in np.arange(nelements):
        f.write("{:3d}  {:2s}  {:5.2f}  {:10s}  {:12.8f}\n".
                format(index[i], symbol[i], dex[i], name[i], mass[i]))


def readatomic(afile):
  """
  Read an elemental (atomic) composition file.

  Parameters
  ----------
  afile: String
    File with atomic composition.

  Returns
  -------
  anum: 1D integer ndarray
     Atomic number (except for Deuterium, which has anum=0).
  symbol: 1D string ndarray
     Elemental chemical symbol.
  dex: 1D float ndarray
     Logarithmic number-abundance, scaled to log(H) = 12.
  name: 1D string ndarray
     Element names.
  mass: 1D float ndarray
     Elemental mass in amu.

  Uncredited developers
  ---------------------
  Jasmina Blecic
  """
  # Allocate arrays:
  nelements = 84  # Fixed number
  anum   = np.zeros(nelements, np.int)
  symbol = np.zeros(nelements, '|U2')
  dex    = np.zeros(nelements, np.double)
  name   = np.zeros(nelements, '|U20')
  mass   = np.zeros(nelements, np.double)

  # Open-read file:
  with open(afile, 'r') as f:
      # Read-discard first two lines (header):
      f.readline()
      f.readline()
      # Store data into the arrays:
      for i in np.arange(nelements):
          anum[i], symbol[i], dex[i], name[i], mass[i] = f.readline().split()

  return anum, symbol, dex, name, mass


def makepreatm(pressure, temp, afile, elements, species, patm):
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
  index, symbol, dex, name, mass = readatomic(afile)
  # Take only the elements we need:
  iatoms = np.in1d(symbol, elements)

  # Lisf of fractional number of molecules relative to hydrogen:
  nfrac = (10.0**(dex-dex[np.where(symbol=="H")]))[iatoms].tolist()

  # Write pre-atm file:
  f = open(patm, 'w')

  # Pre-atm header with basic instructions
  f.write("# TEA pre-atmosphere file.\n"
          "# Units: pressure (bar), temperature (K), abundance (unitless).\n\n")

  # Load defaults:
  sdefaults = {}
  with open(pc.ROOT + "inputs/TEA_gdata_defaults.txt", "r") as d:
      for line in d:
          sdefaults[line.split()[0]] = line.split()[1]

  # Check species names:
  for i in np.arange(len(species)):
      if species[i].find("_") < 0:
          species[i] = sdefaults[species[i]]

  # Write species names
  f.write('#SPECIES\n{:s}\n\n'.format(" ".join(species)))

  # Write TEA data:
  f.write("#TEADATA\n#Pressure          Temp  " +
          "  ".join(["{:>12s}".format(atom) for atom in symbol[iatoms]]) + "\n")
  for i in np.arange(nlayers):
      f.write("{:10.4e}     {:>8.2f}  ".format(pressure[i], temp[i]))
      f.write("  ".join(["{:12.6e}".format(abun) for abun in nfrac]) + "\n")
  f.close()


def abundances(atmfile, pressure, temperature, species, elements=None,
               quniform=None, punits="bar", xsolar=1.0,
               solar=pc.ROOT+"inputs/AsplundEtal2009.txt",
               log=None, verb=0):
  """
  Wrapper to compute atmospheric abundaces for given pressure and
  temperature profiles with either uniform abundances or TEA.

  Parameters
  ----------
  atmfile: String
     Output file where to save the atmospheric model.
  pressure: 1D float ndarray
     Atmospheric pressure profile (barye).
  temperature: 1D float ndarray
     Atmospheric temperature profile (Kelvin).
  species: 1D string list
     Atmospheric composition.
  elements: 1D strings list
     Elemental composition.
  quniform: 1D float ndarray
     If not None, the species abundances at all layers.
  punits: String
     Output pressure units.
  xsolar: Float
     Metallicity enhancement factor.
  solar: String
     Solar elemental abundances file.
  log: Log object
     Screen-output log handler.
  verb: Integer
     Verbosity level.

  Returns
  -------
  q: 2D float ndarray
     Atmospheric abundances (volume mixing fraction).

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
  >>> Q = pa.abundances("pbtea.atm", press, temp, species)
  """
  if log is None:
      log = mu.Log(verb=verb)
  # Uniform-abundances profile:
  if quniform is not None:
      q = uniform(pressure, temperature, species, quniform, punits, log,
                  atmfile)
      log.head("\nProduced uniform-abundances atmospheric file: '{:s}'.".
               format(atmfile))
      return q

  # TEA abundances:
  log.head("\nRun TEA to compute thermochemical-equilibrium abundances.")
  # Prep up files:
  atomicfile, patm = "PBatomicfile.tea", "PBpreatm.tea"
  makeatomic(solar, atomicfile, xsolar, swap=None)
  if elements is None:
     elements, dummy = stoich(species)
  specs = elements + list(np.setdiff1d(species, elements))
  makepreatm(pressure/pt.u(punits), temperature, atomicfile, elements,
             specs, patm)
  # Run TEA:
  pt.make_tea(abun_file=atomicfile)
  proc = subprocess.Popen([pc.ROOT + "modules/TEA/tea/runatm.py", patm, "TEA"])
  proc.communicate()
  # Reformat the TEA output into the pyrat format:
  TEA2pyrat("TEA.tea", atmfile, species)
  os.remove(atomicfile)
  os.remove(patm)
  os.remove("TEA.cfg")
  os.remove("TEA.tea")
  log.head("Produced TEA atmospheric file '{:s}'.".format(atmfile))
  q = io.read_atm(atmfile)[4]
  return q


def TEA2pyrat(teafile, atmfile, req_species):
  """
  Format a TEA atmospheric file into a Pyrat atmospheric file.

  Paramters
  ---------
  teafile:  String
     Input TEA atmospheric file.
  atmfile:  String
     Output Pyrat atmospheric file.
  req_species: List of strings
     The requested species for output.
  """
  # Open read the TEA file:
  f = open(teafile, "r")
  tea = np.asarray(f.readlines())
  f.close()

  # Line-indices where the species and data are:
  ispec = np.where(tea == "#SPECIES\n")[0][0] + 1  # Species list
  idata = np.where(tea == "#TEADATA\n")[0][0] + 2  # data starting line

  # Read and clean species names:
  species = tea[ispec].split()
  nspecies    = len(species)
  nreqspecies = len(req_species)
  nreqspecies = len(req_species)
  for i in range(nspecies):
      species[i] = species[i].rpartition("_")[0]

  # Species indices corresponding to req_species:
  sindex = np.zeros(len(req_species), int)
  for i in range(len(req_species)):
      sindex[i] = species.index(req_species[i])

  # Extract per-layer data:
  nlayers = len(tea) - idata
  temperature = np.zeros(nlayers)
  pressure    = np.zeros(nlayers)
  abundance   = np.zeros((nlayers, nreqspecies))
  for i in range(nlayers):
      data = np.asarray(tea[idata+i].split(), np.double)
      pressure[i], temperature[i] = data[0:2]
      abundance[i,:] = data[2:][sindex]

  # TEA pressure units are always bars:
  punits = "bar"
  pressure = pressure * pt.u(punits)  # Set in barye units
  # File header:
  header = "# TEA atmospheric file formatted for Pyrat.\n\n"

  io.write_atm(atmfile, pressure, temperature, req_species, abundance,
               punits, header)


def pressure(ptop, pbottom, nlayers, units="bar", log=None, verb=0):
  """
  Compute a log-scale pressure profile.

  Parameters
  ----------
  ptop: String or Float
     Pressure at the top of the atmosphere. If string, may contain the units.
  pbottom: String or Float
     Pressure at the bottom of the atmosphere. If string, may contain the units.
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
  ptop    = pt.get_param('ptop',    ptop,    units, log, gt=0.0)
  pbottom = pt.get_param('pbottom', pbottom, units, log, gt=0.0)
  if ptop >= pbottom:
      log.error('Bottom-layer pressure ({:.2e} {:s}) must be higher than the '
                'top-layer pressure ({:.2e} {:s}).'.
                format(pbottom/pt.u(units), units, ptop/pt.u(units), units))

  # Create pressure array in barye (CGS) units:
  press = np.logspace(np.log10(ptop), np.log10(pbottom), nlayers)
  log.head("Creating {:d}-layer atmospheric model between {:.1e} and "
      "{:.1e} bar.".format(nlayers, ptop/pt.u(units), pbottom/pt.u(units)))
  return press


def temperature(tmodel, pressure=None, rstar=None, tstar=None, tint=100.0,
                gplanet=None, smaxis=None, runits="cm", nlayers=None,
                log=None, tparams=None):
  """
  Temperature profile wrapper.

  Parameters
  ----------
  tmodel: String
      Name of the temperature model.
  pressure: 1D float ndarray
      Atmospheric pressure profile in barye units.
  rstar: String or float
      Stellar radius. If string, may contain the units.
  tstar: String or float
      Stellar temperature in kelvin degrees.
  tint: String or float
      Planetary internal temperature in kelvin degrees.
  gplanet: String or float
      Planetary atmospheric temperature in cm s-2.
  smaxis: String or float
      Orbital semi-major axis. If string, may contain the units.
  runits: String
      Default units for rstar and smaxis.  Available units are: A, nm, um,
      mm, cm (default), m, km, au, pc, rearth, rjup, rsun.
  nlayers: Integer
      Number of pressure layers.
  log: Log object
      Screen-output log handler.
  tparams: 1D float ndarray
      Temperature model parameters. If None, return a tuple with the
      temperature model, its arguments, and the number or required parameters.

  Returns
  -------
  If tparams is not None:
      temperature: 1D float ndarray
          The evaluated atmospheric temperature profile.
  If tparams is None:
      temp_model: Callable
          The atmospheric temperature model.
      targs: List
          The list of additional arguments (besides the model parameters).
      ntpars: Integer
          The expected number of model parameters.

  Examples
  --------
  >>> import pyratbay.atmosphere as pa
  >>> nlayers = 11
  >>> # Isothermal profile:
  >>> temp_iso = pa.temperature("isothermal", tparams=1500.0, nlayers=nlayers)
  >>> print(temp_iso)
  [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]
  >>> # Three-channel Eddington-approximation profile:
  >>> pressure = pa.pressure(1e-8, 1e2, nlayers, "bar")
  >>> tparams = np.array([-1.5, -0.8, -0.8, 0.5, 1.0])
  >>> temp = pa.temperature('tcea', pressure, rstar="0.756 rsun", tstar=5040,
  >>>     tint=100.0, gplanet=2200.0, smaxis="0.031 au", tparams=tparams)
  >>> print(temp)
  [1047.04157312 1047.04189805 1047.04531644 1047.08118784 1047.45648563
   1051.34469989 1088.69956369 1311.86379107 1640.12857767 1660.02396061
   1665.30121021]
  """
  if log is None:
      log = mu.Log()

  if tmodel == 'tcea':
      # Parse inputs:
      rstar   = pt.get_param('rstar',   rstar,   runits,   log, gt=0.0)
      tstar   = pt.get_param('tstar',   tstar,   'kelvin', log, gt=0.0)
      tint    = pt.get_param('tint',    tint,    'kelvin', log, gt=0.0)
      gplanet = pt.get_param('gplanet', gplanet, 'none',   log, gt=0.0)
      smaxis  = pt.get_param('smaxis',  smaxis,  runits,   log, gt=0.0)
      # Define model and arguments:
      targs  = [pressure, rstar, tstar, tint, gplanet, smaxis]
      temp_model = tmodels.TCEA(*targs)
      ntpars = 5
  elif tmodel == 'isothermal':
      targs  = [nlayers]
      temp_model = tmodels.Isothermal(*targs)
      ntpars = 1
  elif tmodel == 'madhu':
      targs = [pressure]
      temp_model = tmodels.Madhu(*targs)
      ntpars = 6
  else:
      log.error("Invalid input temperature model '{:s}'.  Select from: {}".
          format(tmodel, pc.tmodels))

  if tparams is None:
      return temp_model, targs, ntpars
  else:
      temperature = temp_model(tparams)
      log.head('\nComputed {:s} temperature model.'.format(tmodel))
      return temperature


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
  radius = si.cumtrapz((-pc.k*sc.N_A * temperature) / (mu*g), np.log(pressure))
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


def hydro_m(pressure, temperature, mu, M, p0, r0):
  """
  Calculate radii using the hydrostatic-equilibrium equation considering
  a variable gravity: g=G*M/r**2

  Parameters
  ----------
  pressure: 1D float ndarray
     Atmospheric pressure for each layer (in barye).
  temperature: 1D float ndarray
     Atmospheric temperature for each layer (in K).
  mu: 1D float ndarray
     Mean molecular mass for each layer (in g mol-1).
  M: Float
     Planetary mass (in g).
  p0: Float
     Reference pressure level (in barye) where radius(p0) = r0.
  r0: Float
     Reference radius level (in cm) corresponding to p0.

  Returns
  -------
  radius: 1D float ndarray
     Radius for each layer (in cm).

  Examples
  --------
  >>> import pyratbay.atmosphere as pa
  >>> import pyratbay.constants  as pc
  >>> nlayers = 11
  >>> pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
  >>> temperature = pa.tmodels.Isothermal(nlayers)(1500.0)
  >>> mu = np.tile(2.3, nlayers)
  >>> Mp = 1.0 * pc.mjup
  >>> r0 = 1.0 * pc.rjup
  >>> p0 = 1.0 * pc.bar
  >>> # Radius profile in Jupiter radii:
  >>> radius = pa.hydro_m(pressure, temperature, mu, Mp, p0, r0) / pc.rjup
  >>> print(radius)
  [1.05973436 1.05188019 1.04414158 1.036516   1.029001   1.02159419
   1.01429324 1.00709591 1.         0.99300339 0.986104  ]
  """
  # Apply the HE equation:
  I = si.cumtrapz((pc.k*sc.N_A*temperature)/(pc.G*mu*M), np.log(pressure))
  I = np.concatenate(([0.0], I))

  # Find current radius at p0:
  radinterp = sip.interp1d(pressure, I, kind='slinear')
  I0 = radinterp(p0)
  # Set: radius(p0) = r0
  radius = 1.0/(I - I0 + 1/r0)

  return radius


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
          # New element:
          if char.isupper():
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


def readmol(file):
  """
  Read a molecules file to extract their ID, symbol, mass, and diameter.

  Parameters
  ----------
  file: String
     The molecule file path.

  Returns
  -------
  molID: 1D integer ndarray
     The molecules' ID.
  symbol: 1D string ndarray
     The molecule's name.
  mass: 1D float ndarray
     The mass of the molecules (in g mol-1).
  diam: 1D float ndarray
     The collisional diameter of the molecules (in Angstrom).

  Notes
  -----
  In all truthfulness, these are species, not only molecules, as the
  file also contain elemental particles.
  """
  with open(file, 'r') as molfile:
      # Skip comment and blank lines:
      line = molfile.readline().strip()
      while line == '' or line.startswith('#'):
          line = molfile.readline().strip()

      # Allocate outputs:
      molID  = [] # Molecule ID
      symbol = [] # Molecule symbol
      mass   = [] # Molecule mass
      diam   = [] # Molecule diameter

      # Read in values:
      while line != '' and not line.startswith('#'):  # Start reading species
          molinfo = line.split()
          # Extract info:
          molID .append(  int(molinfo[0]))
          symbol.append(      molinfo[1] )
          mass  .append(float(molinfo[2]))
          diam  .append(float(molinfo[3]))
          line = molfile.readline().strip()  # Read next line

  molID  = np.asarray(molID)
  symbol = np.asarray(symbol)
  mass   = np.asarray(mass)
  diam   = np.asarray(diam)

  return molID, symbol, mass, diam


def meanweight(abundances, species, molfile=pc.ROOT+'inputs/molecules.dat'):
  """
  Calculate the mean molecular weight (a.k.a. mean molecular mass)
  for the given abundances composition.

  Parameters
  ----------
  abundances: 2D float iterable
      Species mol-mixing-fraction array of shape [nlayers,nmol].
  species: 1D string iterable
      Species names.
  molfile: String
      A molecules file with the species info.

  Returns
  -------
  mu: 1D float ndarray
      Mean molecular weight at each layer for the input abundances.

  Examples
  --------
  >>> import pyratbay.atmosphere as pa
  >>> species     = ['H2', 'He', 'H2O', 'CO', 'CO2', 'CH4']
  >>> abundances  = [[0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]]
  >>> mu = pa.meanweight(abundances, species)
  >>> print(mu)
  [2.31928918]
  """
  molID, symbol, mass, diam = readmol(molfile)
  molmass = np.array([mass[symbol==spec][0] for spec in species])
  return np.sum(np.atleast_2d(abundances)*molmass, axis=1)


def IGLdensity(abundances, pressure, temperature):
  """
  Use the Ideal gas law to calculate the density in molecules cm-3.

  Parameters
  ----------
  abundances: 2D float ndarray
     Species volume mixing ratio profiles, of shape [nlayers,nmol].
  pressure: 1D ndarray
     Atmospheric pressure profile (in barye units).
  temperature: 1D ndarray
     Atmospheric temperature (in kelvin).

  Returns
  -------
  density: 2D float ndarray
     Atmospheric density in molecules per centimeter^3.

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
  >>> dens = pa.IGLdensity(qprofiles, pressure, temperature)
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

