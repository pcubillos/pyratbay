# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np
import scipy.interpolate as sip

from .. import tools      as pt
from .. import constants  as pc
from .. import atmosphere as pa


def readatm(pyrat):
  """
  Read the atmospheric file, store variables in pyrat.
  """
  # Check atmfile:
  pyrat.log.msg("\nReading atmospheric file: '{:s}'.".format(pyrat.atmfile))

  # User-input atmospheric-data object:
  atm = pyrat.inputs.atm

  try:
      atm_inputs = pa.readatm(pyrat.atmfile)
  except ValueError as e:
      # TBD: handle properly the error message
      pyrat.log.error("Atmosphere file has an unexpected line: \n'{:s}'".
                      format(line))

  # Atmospheric-file units, species, and profiles:
  punits, tunits, qunits, runits = atm_inputs[0]
  pyrat.mol.name = atm_inputs[1]
  atm.press, atm.temp, atm.q, atm.radius = atm_inputs[2:]

  # Warnings:
  atm.punits = pt.defaultp(punits, 'barye', "Undefined pressure units in "
      "the input atmospheric file.  Assumed to be '{:s}'.", pyrat.log)
  atm.tunits = pt.defaultp(tunits, 'kelvin', "Undefined temperature units in "
      "the input atmospheric file.  Assumed to be '{:s}'.", pyrat.log)
  if atm.radius is not None:
      atm.runits = pt.defaultp(runits, 'cm', "Undefined radius units in the "
          "input atmospheric file.  Assumed to be '{:s}'.", pyrat.log)
  atm.qunits = pt.defaultp(qunits, 'number', "Undefined abundance units in "
      "the input atmospheric file.  Assumed to be '{:s}'.", pyrat.log)

  atm.nlayers, pyrat.mol.nmol = np.shape(atm.q)
  pyrat.log.msg("Species list: \n  {:s}".format(str(pyrat.mol.name)),
                verb=2, indent=2, si=4)

  txt = "volume" if (qunits == "number") else "mass"
  pyrat.log.msg("Abundances are given by {:s} ({:s} mixing ratio).".
      format(atm.qunits, txt), verb=2, indent=2)
  pyrat.log.msg("Unit factors: radius: {}, pressure: {:s}, temperature: {:s}".
      format(atm.runits, atm.punits, atm.tunits), verb=2, indent=2)

  # Read molecular constant values:
  getconstants(pyrat)

  # Store values in CGS system of units:
  if atm.radius is not None:
      atm.radius *= pt.u(atm.runits)
  atm.press *= pt.u(atm.punits)
  atm.temp  *= pt.u(atm.tunits)

  pyrat.log.msg("Number of layers in the input atmospheric file: {:d}".
                format(atm.nlayers), verb=2, indent=2)
  pyrat.log.msg("Atmospheric file pressure limits: {:.2e}--{:.2e} {:s}.".
      format(atm.press[ 0]/pt.u(atm.punits),
             atm.press[-1]/pt.u(atm.punits), atm.punits), verb=2, indent=2)

  # Calculate the mean molecular mass per layer:
  atm.mm = np.sum(atm.q*pyrat.mol.mass, axis=1)
  pyrat.log.msg("Typical mean molecular mass: {:.3f} g mol-1.".
                format(np.median(atm.mm)), verb=2, indent=2)

  # Store the abundance as volume mixing ratio:
  if atm.qunits == "mass":
      atm.q = atm.q * atm.mm / pyrat.mol.mass

  # Calculate number density profiles for each molecule (in molecules cm-3):
  atm.d = pa.IGLdensity(atm.q, atm.press, atm.temp)

  pyrat.log.msg("Read atmosphere done.")


def getconstants(pyrat):
  """
  Set molecules constant values (ID, mass, radius).
  """
  # Read file with molecular info:
  pyrat.log.msg("Taking species constant parameters from: '{:s}'.".
                format(pyrat.molfile), verb=2, indent=2)
  molID, symbol, mass, diam = pa.readmol(pyrat.molfile)

  # Check that all atmospheric species are listed in molfile:
  absent = np.setdiff1d(pyrat.mol.name, symbol)
  if len(absent) > 0:
    pyrat.log.error("These species: {:s} are not listed in the molecules "
                    "info file: {:s}.\n".format(str(absent), pyrat.molfile))

  # Set molecule's values:
  pyrat.mol.ID     = np.zeros(pyrat.mol.nmol, np.int)
  pyrat.mol.symbol = np.zeros(pyrat.mol.nmol, "U15")
  pyrat.mol.mass   = np.zeros(pyrat.mol.nmol)
  pyrat.mol.radius = np.zeros(pyrat.mol.nmol)

  pyrat.log.msg("Molecule   ID   Radius  Mass\n"
                "                (A)     (gr/mol)", verb=2, indent=4)
  for i in np.arange(pyrat.mol.nmol):
    # Find the molecule in the list:
    imol = np.where(symbol == pyrat.mol.name[i])[0]
    # Set molecule ID:
    pyrat.mol.ID[i]     = molID [imol]
    # Set molecule symbol:
    pyrat.mol.symbol[i] = symbol[imol][0]
    # Set molecule mass:
    pyrat.mol.mass[i]   = mass  [imol]
    # Set molecule collision radius:
    pyrat.mol.radius[i] = 0.5 * diam[imol] * pc.A
    pyrat.log.msg("{:>10s}:  {:3d}  {:.3f}  {:8.4f}".format(pyrat.mol.name[i],
                  pyrat.mol.ID[i], pyrat.mol.radius[i]/pc.A,
                  pyrat.mol.mass[i]), verb=2, indent=2)


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
      temp = pyrat.atm.tmodel(pyrat.atm.tpars, *pyrat.atm.targs)
  else:
      temp = pyrat.atm.temp

  # Check that the dimensions match:
  if np.size(temp) != np.size(pyrat.atm.temp):
      pyrat.log.error("The temperature array size ({:d}) doesn't match the "
                      "Pyrat's temperature size ({:d}).".
                      format(np.size(temp), np.size(pyrat.atm.temp)))

  # Check temperature boundaries:
  errorlog = ("One or more input temperature values lies out of the {:s} "
      "temperature boundaries (K): [{:6.1f}, {:6.1f}].\nHalted calculation.")
  if pyrat.ex.extfile is not None:
      if np.any(temp > pyrat.ex.tmax) or np.any(temp < pyrat.ex.tmin):
          pyrat.log.warning(errorlog.format("tabulated EC", pyrat.ex.tmin,
                                                            pyrat.ex.tmax))
          return 0
  else:
      if (pyrat.lt.ntransitions > 0 and
         (np.any(temp > pyrat.lt.tmax) or np.any(temp < pyrat.lt.tmin))):
          pyrat.log.warning(errorlog.format("line-transition", pyrat.lt.tmin,
                                                               pyrat.lt.tmax))
          return 0
      if (pyrat.cs.nfiles > 0 and
         (np.any(temp > pyrat.cs.tmax) or np.any(temp < pyrat.cs.tmin))):
          pyrat.log.warning(errorlog.format("cross-section", pyrat.cs.tmin,
                                                             pyrat.cs.tmax))
          return 0

  # Recompute abundance profiles:
  q0 = np.copy(pyrat.atm.qbase)
  if abund is not None:
      pass
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
  else:
      # Compute gplanet from mass and radius if necessary/possible:
      if pyrat.phy.gplanet is None and \
         pyrat.phy.mplanet is not None and \
         pyrat.phy.rplanet is not None:
          pyrat.phy.gplanet = pc.G * pyrat.phy.mplanet / pyrat.phy.rplanet**2

      # Check that the gravity variable is exists:
      if pyrat.phy.rplanet is None:
          pyrat.log.error("Undefined reference planetary radius (rplanet). "
           "Either provide the radius profile for the layers or the rplanet.")
      if pyrat.phy.gplanet is None:
          pyrat.log.error("Undefined atmospheric gravity (gplanet). "
           "Either provide the radius profile for the layers, the surface "
           "gravity, or the planetary mass (mplanet).")
      if pyrat.refpressure is None:
          pyrat.log.error("Undefined reference pressure level (refpressure). "
           "Either provide the radius profile for the layers or refpressure.")
      pyrat.atm.radius = pyrat.hydro(pyrat.atm.press, pyrat.atm.temp,
                            pyrat.atm.mm, pyrat.phy.gplanet, pyrat.phy.mplanet,
                            pyrat.refpressure, pyrat.phy.rplanet)
  # Check radii lie within Hill radius:
  rtop = np.where(pyrat.atm.radius > pyrat.phy.rhill)[0]
  if np.size(rtop) > 0:
    pyrat.atm.rtop = rtop[-1] + 1
    pyrat.log.warning("The atmospheric pressure array extends beyond the "
       "Hill radius ({:.3e} km) at pressure {:.3e} bar (layer #{:d}).  "
       "Extinction beyond this layer will be neglected.".
        format(pyrat.phy.rhill/pc.km, pyrat.atm.press[pyrat.atm.rtop]/pc.bar,
               pyrat.atm.rtop))
  else:
    pyrat.atm.rtop = 0

  # Partition function:
  for db in pyrat.lt.db:            # For each Database
      for j in np.arange(db.niso):  # For each isotope in DB
          zinterp = sip.interp1d(db.temp, db.z[j], kind='slinear')
          pyrat.iso.z[db.iiso+j] = zinterp(pyrat.atm.temp)
