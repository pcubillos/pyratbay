# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np
import scipy.interpolate as sip

from .. import tools      as pt
from .. import constants  as pc
from .. import atmosphere as pa
from .. import io as io


def read_atm(pyrat):
  """
  Read the atmospheric file, store variables in pyrat.
  """
  # Check atmfile:
  pyrat.log.head("\nReading atmospheric file: '{:s}'.".
      format(pyrat.atm.atmfile))

  # User-input atmospheric-data object:
  atm    = pyrat.atm
  atm_in = pyrat.inputs.atm

  with pt.log_error(pyrat.log):
      atm_inputs = io.read_atm(pyrat.atm.atmfile)

  # Atmospheric-file units, species, and profiles:
  punits, tunits, qunits, runits = atm_inputs[0]
  pyrat.mol.name = atm_inputs[1]
  pyrat.mol.nmol = len(pyrat.mol.name)

  # Read molecular constant values:
  get_constants(pyrat)

  # Store values in CGS system of units:
  atm_in.press  = atm_inputs[2] * pt.u(punits)
  atm_in.temp   = atm_inputs[3] * pt.u(tunits)
  atm_in.q      = atm_inputs[4]
  if atm_inputs[5] is not None:
      atm_in.radius = atm_inputs[5] * pt.u(runits)
  atm_in.nlayers = len(atm_in.press)

  # Calculate the mean molecular mass per layer:
  if qunits == 'mass':
      atm_in.mm = 1.0/np.sum(atm_in.q/pyrat.mol.mass, axis=1)
  elif qunits == 'number':
      atm_in.mm = np.sum(atm_in.q*pyrat.mol.mass, axis=1)
  else:
      pyrat.log.error("Invalid input abundance units ('{:s}').".format(qunits))

  # Store the abundance as volume mixing ratio:
  if qunits == "mass":
      atm_in.q = atm_in.q * atm_in.mm / pyrat.mol.mass

  # Calculate number density profiles for each molecule (in molecules cm-3):
  atm_in.d = pa.IGLdensity(atm_in.q, atm_in.press, atm_in.temp)

  pyrat.log.msg("Species list:\n  {:s}".format(str(pyrat.mol.name)),
                indent=2, si=4)

  txt = "volume" if (qunits == "number") else "mass"
  pyrat.log.msg("Abundances are given by {:s} ({:s} mixing ratio).".
      format(qunits, txt), indent=2)
  pyrat.log.msg("Unit factors: radius: {}, pressure: {:s}, temperature: {:s}".
      format(runits, punits, tunits), indent=2)

  pyrat.log.msg("Number of layers in the input atmospheric file: {:d}".
      format(atm_in.nlayers), indent=2)
  pyrat.log.msg("Atmospheric file pressure limits: {:.2e}--{:.2e} {:s}.".
      format(atm_in.press[ 0]/pt.u(atm.punits),
             atm_in.press[-1]/pt.u(atm.punits), atm.punits), indent=2)

  pyrat.log.msg("Typical mean molecular mass: {:.3f} g mol-1.".
      format(np.median(atm_in.mm)), indent=2)

  pyrat.log.head("Read atmosphere done.")


def get_constants(pyrat):
  """
  Set molecules constant values (ID, mass, radius).
  """
  # Read file with molecular info:
  pyrat.log.msg("Taking species constant parameters from: '{:s}'.".
                format(pyrat.mol.molfile), indent=2)
  molID, symbol, mass, diam = pa.readmol(pyrat.mol.molfile)

  # Check that all atmospheric species are listed in molfile:
  absent = np.setdiff1d(pyrat.mol.name, symbol)
  if len(absent) > 0:
      pyrat.log.error("These species: {:s} are not listed in the molecules "
          "info file: {:s}.\n".format(str(absent), pyrat.mol.molfile))

  # Set molecule's values:
  pyrat.mol.ID     = np.zeros(pyrat.mol.nmol, np.int)
  pyrat.mol.symbol = np.zeros(pyrat.mol.nmol, 'U15')
  pyrat.mol.mass   = np.zeros(pyrat.mol.nmol)
  pyrat.mol.radius = np.zeros(pyrat.mol.nmol)

  pyrat.log.msg('Molecule   ID   Radius  Mass\n'
                '                (A)     (gr/mol)', indent=4)
  for i in range(pyrat.mol.nmol):
      # Find the molecule in the list:
      imol = np.where(symbol == pyrat.mol.name[i])[0]
      # Set molecule ID, name, mass, and collision radius:
      pyrat.mol.ID[i]     = molID [imol]
      pyrat.mol.symbol[i] = symbol[imol][0]
      pyrat.mol.mass[i]   = mass  [imol]
      pyrat.mol.radius[i] = 0.5*diam[imol] * pc.A
      pyrat.log.msg("{:>10s}:  {:3d}  {:.3f}  {:8.4f}".format(pyrat.mol.name[i],
                    pyrat.mol.ID[i], pyrat.mol.radius[i]/pc.A,
                    pyrat.mol.mass[i]), indent=2)


def reloadatm(pyrat, temp=None, abund=None, radius=None):
  """
  Parameters
  ----------
  pyrat: A Pyrat instance
  temp: 1D float ndarray
      Layer's temperature profile (in Kelvin) sorted from top to bottom.
  abund: 2D float ndarray
      Species mole mixing ratio profiles [nlayers, nmol].
  radius: 1D float ndarray
      Layer's altitude profile (in cm), same order as temp.
  """
  # Recompute temperature profile:
  if temp is not None:
      # Need to null tpars since it does not represent temp anymore
      pyrat.atm.tpars = None
  elif pyrat.atm.tpars is not None:
      temp = pyrat.atm.tmodel(pyrat.atm.tpars)
  else:
      temp = pyrat.atm.temp

  # Check that the dimensions match:
  if np.size(temp) != np.size(pyrat.atm.temp):
      pyrat.log.error("The temperature array size ({:d}) doesn't match the "
                      "Pyrat's temperature size ({:d}).".
                      format(np.size(temp), np.size(pyrat.atm.temp)))

  # Check temperature boundaries:
  errorlog = ("One or more input temperature values lies out of the {:s} "
              "temperature boundaries (K): [{:6.1f}, {:6.1f}].")
  if pyrat.ex.extfile is not None:
      if np.any(temp > pyrat.ex.tmax) or np.any(temp < pyrat.ex.tmin):
          pyrat.log.warning(errorlog.format('tabulated extinction-coefficient',
                                            pyrat.ex.tmin, pyrat.ex.tmax))
          return 0
  elif pyrat.lt.ntransitions > 0:
      if np.any(temp > pyrat.lt.tmax) or np.any(temp < pyrat.lt.tmin):
          pyrat.log.warning(errorlog.format('line-transition',
                                            pyrat.lt.tmin, pyrat.lt.tmax))
          return 0
  if pyrat.cs.nfiles > 0:
      if np.any(temp > pyrat.cs.tmax) or np.any(temp < pyrat.cs.tmin):
          pyrat.log.warning(errorlog.format('cross-section',
                                            pyrat.cs.tmin, pyrat.cs.tmax))
          return 0

  # Recompute abundance profiles:
  q0 = np.copy(pyrat.atm.qbase)
  if abund is not None:
      pyrat.atm.molpars = None
  elif pyrat.atm.molpars is not None:
      abund = pa.qscale(q0, pyrat.mol.name, pyrat.atm.molmodel,
          pyrat.atm.molfree, pyrat.atm.molpars, pyrat.atm.bulk,
          iscale=pyrat.atm.ifree, ibulk=pyrat.atm.ibulk,
          bratio=pyrat.atm.bulkratio, invsrat=pyrat.atm.invsrat)
  else:
      abund = q0
  if np.shape(abund) != np.shape(pyrat.atm.q):
      pyrat.log.error("The shape of the abundances array {:s} doesn't match "
                      "the shape of the Pyrat's abundance size {:s}".
                      format(str(np.shape(abund)), str(np.shape(pyrat.atm.q))))

  # Update values:
  pyrat.atm.temp = temp
  pyrat.atm.q    = abund

  # Mean molecular mass:
  pyrat.atm.mm = np.sum(pyrat.atm.q*pyrat.mol.mass, axis=1)

  # Number density (molecules cm-3):
  pyrat.atm.d = pa.IGLdensity(pyrat.atm.q, pyrat.atm.press, pyrat.atm.temp)

  # Take radius if provided, else use hydrostatic-equilibrium equation:
  if radius is not None:
      pyrat.atm.radius = radius
  elif (pyrat.inputs.atm.radius is not None
      and (pyrat.atm.refpressure is None or pyrat.phy.rplanet is None)):
          pyrat.atm.radius = pyrat.inputs.atm.radius
  else:
      # Compute gplanet from mass and radius if necessary/possible:
      if pyrat.phy.gplanet is None and     \
         pyrat.phy.mplanet is not None and \
         pyrat.phy.rplanet is not None:
          pyrat.phy.gplanet = pc.G * pyrat.phy.mplanet / pyrat.phy.rplanet**2

      # Check that the gravity variable is exists:
      if pyrat.phy.rplanet is None:
          pyrat.log.error('Undefined reference planetary radius (rplanet). '
           'Either provide the radius profile or define rplanet.')
      if pyrat.phy.gplanet is None:
          pyrat.log.error('Undefined atmospheric gravity (gplanet). Either '
           'provide the radius profile, define gplanet, or define mplanet.')
      if pyrat.atm.refpressure is None:
          pyrat.log.error('Undefined reference pressure level (refpressure). '
           'Either provide the radius profile or define refpressure.')
      pyrat.atm.radius = pyrat.hydro(pyrat.atm.press, pyrat.atm.temp,
                            pyrat.atm.mm, pyrat.phy.gplanet, pyrat.phy.mplanet,
                            pyrat.atm.refpressure, pyrat.phy.rplanet)

  # Check radii lie within Hill radius:
  pyrat.atm.rtop = pt.ifirst(pyrat.atm.radius < pyrat.phy.rhill, default_ret=0)
  if pyrat.atm.rtop > 0:
      pyrat.log.warning("The atmospheric pressure array extends beyond the "
          "Hill radius ({:.3e} km) at pressure {:.3e} bar (layer {:d}).  "
          "Extinction beyond this layer will be neglected.".
          format(pyrat.phy.rhill/pc.km, pyrat.atm.press[pyrat.atm.rtop]/pc.bar,
                 pyrat.atm.rtop))

  # Partition function:
  for db in pyrat.lt.db:            # For each Database
      for j in np.arange(db.niso):  # For each isotope in DB
          zinterp = sip.interp1d(db.temp, db.z[j], kind='slinear')
          pyrat.iso.z[db.iiso+j] = zinterp(pyrat.atm.temp)
