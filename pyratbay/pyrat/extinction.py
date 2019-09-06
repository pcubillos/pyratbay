# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import ctypes
import multiprocessing as mpr

import numpy as np
import scipy.interpolate as sip

from .. import tools     as pt
from .. import constants as pc
from .. import io        as io

sys.path.append(pc.ROOT + "pyratbay/lib/")
import extcoeff as ec


def exttable(pyrat):
  """
  Handle extinction-coefficient table (read/calculate/write file).
  """
  # ID list of species with isotopes:
  if np.size(np.where(pyrat.iso.imol>=0)[0]) == 0:
      pyrat.ex.nmol = 0
  else:
      imols = pyrat.iso.imol[np.where(pyrat.iso.imol>=0)]
      pyrat.ex.molID = pyrat.mol.ID[np.unique(imols)]
      pyrat.ex.nmol  = len(pyrat.ex.molID)

  # If the extinction file was not defined, skip this step:
  if pyrat.ex.extfile is None:
      pyrat.log.head("\nNo extinction-coefficient table requested.")
      return

  # If the extinction file exists, read it:
  if pt.isfile(pyrat.ex.extfile) == 1:
      pyrat.log.head("\nReading extinction-coefficient table file(s):")
      for extfile in pyrat.ex.extfile:
          pyrat.log.head("  '{}'.".format(extfile))
      read_extinction(pyrat)
  # If the extinction file doesn't exist, calculate it:
  else:
      pyrat.log.head("\nGenerating new extinction-coefficient table file:"
                     "\n  '{:s}'.".format(pyrat.ex.extfile[0]))
      calc_extinction(pyrat)


def read_extinction(pyrat):
  """
  Read extinction-coefficient tables from files.
  """
  ex = pyrat.ex
  ex.molID, ex.etable = [], []

  for extfile in ex.extfile:
      edata = io.read_opacity(extfile)
      # Array sizes:
      nmol, ntemp, nlayers, nwave = edata[0]
      if (ex.ntemp is not None and
          (ntemp != ex.ntemp or nlayers != ex.nlayers or nwave != ex.nwave)):
          pyrat.log.error("Shape of the extinction-coefficient file '{:s}' "
              "does not match with previous file shapes.".
              format(os.path.basename(extfile)))
      ex.ntemp, ex.nlayers, ex.nwave = ntemp, nlayers, nwave

      # Molecule ID, temperature (K), pressure (barye), and wavenumber (cm-1):
      molID, temp, press, wn = edata[1]
      if ex.temp is not None and np.any(np.abs(1-temp/ex.temp) > 0.01):
          pyrat.log.error("Tabulated temperature array of file '{:s}' "
              "does not match with the previous temperature arrays.".
              format(os.path.basename(extfile)))
      if ex.press is not None and np.any(np.abs(1-press/ex.press) > 0.01):
          pyrat.log.error("Tabulated pressure array of file '{:s}' "
              "does not match with the previous pressure arrays.".
              format(os.path.basename(extfile)))
      if ex.wn is not None and np.any(np.abs(1-wn/ex.wn) > 0.01):
          pyrat.log.error("Tabulated wavelength array of file '{:s}' "
              "does not match with the previous wavelength arrays.".
              format(os.path.basename(extfile)))
      ex.temp, ex.press, ex.wn = temp, press, wn
      ex.molID += list(molID)
      # Extinction-coefficient table (cm2 molecule-1):
      ex.etable.append(edata[2])

  ex.molID = np.array(ex.molID)
  ex.nmol = len(ex.molID)
  ex.etable = np.vstack(ex.etable)
  if np.ndim(ex.etable) == 3:
      ex.etable = np.expand_dims(ex.etable, axis=0)

  pyrat.log.msg("File(s) have {:d} molecules, {:d} temperature samples, "
                "{:d} layers, and {:d} wavenumber samples.".format(ex.nmol,
                ex.ntemp, ex.nlayers, ex.nwave), indent=2)

  # Set tabulated temperature extrema:
  ex.tmin = np.amin(ex.temp)
  ex.tmax = np.amax(ex.temp)

  with np.printoptions(precision=1):
      str_temp = str(ex.temp)
  with np.printoptions(formatter={'float':'{: .2f}'.format}, threshold=100):
      str_wn   = str(ex.wn)
  with np.printoptions(formatter={'float':'{:.3e}'.format}):
      str_press = str(ex.press/pc.bar)
  pyrat.log.msg(
      "Molecules' IDs: {}\n"
      "Temperatures (K):\n   {}\n"
      "Pressure layers (bar):\n{}\n"
      "Wavenumber array (cm-1):\n   {}".
      format(ex.molID, str_temp, str_press, str_wn), indent=2)

  # Some checks:
  if ex.nwave != pyrat.spec.nwave or np.sum(np.abs(ex.wn-pyrat.spec.wn)) > 0:
      pyrat.log.warning("Wavenumber sampling from extinction-coefficient "
          "table does not match the input wavenumber sampling.  Adopting "
          "tabulated array with {:d} samples, spacing of {:.2f} cm-1, "
          "and ranges [{:.2f}, {:.2f}] cm-1.".
            format(ex.nwave, ex.wn[1]-ex.wn[0], ex.wn[0], ex.wn[-1]))
      # Update wavenumber sampling:
      pyrat.spec.wn     = ex.wn
      pyrat.spec.nwave  = ex.nwave
      pyrat.spec.wnlow  = ex.wn[ 0]
      pyrat.spec.wnhigh = ex.wn[-1]
      pyrat.spec.wnstep = ex.wn[ 1] - ex.wn[0]
      # Keep wavenumber oversampling factor:
      pyrat.spec.ownstep = pyrat.spec.wnstep / pyrat.spec.wnosamp
      pyrat.spec.onwave  = (pyrat.spec.nwave - 1) * pyrat.spec.wnosamp + 1
      pyrat.spec.own     = np.linspace(pyrat.spec.wn[0], pyrat.spec.wn[-1],
                                       pyrat.spec.onwave)

      # Update interpolated stellar spectrum:
      if pyrat.phy.starflux is not None:
          sinterp = sip.interp1d(pyrat.phy.starwn, pyrat.phy.starflux)
          pyrat.spec.starflux = sinterp(pyrat.spec.wn)
      # Update observational variables:
      pyrat.set_filters()


def calc_extinction(pyrat):
  """
  Compute the extinction-coefficient (per species) for a tabulated grid of
  temperatures and pressures, over a wavenumber array.
  """
  # Extinction-coefficient object:
  ex = pyrat.ex

  # Temperature boundaries check:
  if ex.tmin < pyrat.lt.tmin:
      pyrat.log.error('Requested extinction-coefficient table temperature '
          '(tmin={:.1f} K) below the lowest available TLI temperature '
          '({:.1f} K).'.format(ex.tmin, pyrat.lt.tmin))
  if ex.tmax > pyrat.lt.tmax:
      pyrat.log.error('Requested extinction-coefficient table temperature '
          '(tmax={:.1f} K) above the highest available TLI temperature '
          '({:.1f} K).'.format(ex.tmax, pyrat.lt.tmax))

  # Create the temperature array:
  ex.ntemp = int((ex.tmax-ex.tmin)/ex.tstep) + 1
  ex.temp  = np.linspace(ex.tmin, ex.tmin + (ex.ntemp-1)*ex.tstep, ex.ntemp)

  with np.printoptions(formatter={'float':'{:.1f}'.format}):
      pyrat.log.msg("Temperature sample (K):\n {}".format(ex.temp), indent=2)

  # Evaluate the partition function at the given temperatures:
  pyrat.log.msg("Interpolate partition function.", indent=2)
  ex.z = np.zeros((pyrat.iso.niso, ex.ntemp), np.double)
  for i in np.arange(pyrat.lt.ndb):           # For each Database
      for j in np.arange(pyrat.lt.db[i].niso):  # For each isotope in DB
          zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                                 kind='slinear')
          ex.z[pyrat.lt.db[i].iiso+j] = zinterp(ex.temp)

  # Allocate wavenumber, pressure, and isotope arrays:
  ex.wn    = pyrat.spec.wn
  ex.nwave = pyrat.spec.nwave

  ex.press = pyrat.atm.press
  ex.nlayers = pyrat.atm.nlayers

  # Allocate extinction-coefficient array:
  pyrat.log.msg("Calculate extinction coefficient.", indent=2)
  sm_ect = mpr.Array(ctypes.c_double,
                     np.zeros(ex.nmol*ex.ntemp*ex.nlayers*ex.nwave, np.double))
  ex.etable = np.ctypeslib.as_array(sm_ect.get_obj()).reshape(
                             (ex.nmol, ex.ntemp, ex.nlayers, ex.nwave))

  # Multi-processing extinction calculation (in C):
  processes = []
  indices = np.arange(ex.ntemp*ex.nlayers) % pyrat.ncpu  # CPU indices
  for i in np.arange(pyrat.ncpu):
      proc = mpr.Process(target=extinction,           # grid  add
                  args=(pyrat, np.where(indices==i)[0], True, False))
      processes.append(proc)
      proc.start()
  for proc in processes:
      proc.join()

  # Store values in file:
  io.write_opacity(ex.extfile[0], ex.molID, ex.temp, ex.press, ex.wn,
                   pyrat.ex.etable)
  pyrat.log.head("Extinction-coefficient table written to file: '{:s}'.".
                format(ex.extfile[0]), indent=2)


def extinction(pyrat, indices, grid=False, add=False):
  """
  Python multiprocessing wrapper for the extinction-coefficient (EC)
  calculation function for the atmospheric layers or EC grid.

  Parameters
  ----------
  pyrat: Pyrat Object
  indices: 1D integer list
     The indices of the atmospheric layer or EC grid to calculate.
  grid: Bool
     If True, compute EC per species for EC grid.
     If False, compute EC for atmospheric layer.
  add: Bool
     If True, co-add EC contribution (cm-1) from all species.
     If False, keep EC contribution (cm2 molec-1) from each species separated.
  """
  if add:   # Combined or per-species output
      extinct_coeff = np.zeros((1, pyrat.spec.nwave))
  else:
      extinct_coeff = np.zeros((pyrat.ex.nmol, pyrat.spec.nwave))

  if pyrat.iso.iext is None:
      # Get species indices in extinction-coefficient table for the isotopes:
      pyrat.iso.iext = np.zeros(pyrat.iso.niso, np.int)
      for i in np.arange(pyrat.iso.niso):
        if pyrat.iso.imol[i] != -1:
            pyrat.iso.iext[i] = np.where(pyrat.ex.molID ==
                                     pyrat.mol.ID[pyrat.iso.imol[i]])[0][0]
        else:
            pyrat.iso.iext[i] = -1  # Isotope not in atmosphere
            # FINDME: find patch for this case in ec.extinction()

  # Turn off verb of all processes except the first:
  verb = pyrat.log.verb
  pyrat.log.verb = (0 in indices) * verb

  for i,index in enumerate(indices):
      ilayer = index % pyrat.atm.nlayers  # Layer index
      pressure = pyrat.atm.press[ilayer]  # Layer pressure
      molq     = pyrat.atm.q    [ilayer]  # Molecular abundance
      density  = pyrat.atm.d    [ilayer]  # Molecular density
      if grid:  # Take from grid
          itemp = int(index / pyrat.atm.nlayers)  # Temp. index in EC table
          temp   = pyrat.ex.temp[itemp]
          ziso   = pyrat.ex.z[:,itemp]      # Isotopic ratio
          pyrat.log.msg("Extinction-coefficient table: layer {:3d}/{:d}, "
              "iteration {:3d}/{:d}.".format(ilayer+1, pyrat.atm.nlayers, i+1,
                                             len(indices)), indent=2)
      else:     # Take from atmosphere
          temp = pyrat.atm.temp [ilayer]
          ziso = pyrat.iso.z  [:,ilayer]  # Isotopic ratio
          pyrat.log.msg("Calculating extinction at layer {:3d}/{:d} "
              "(T={:6.1f} K, p={:.1e} bar).".format(ilayer+1,
              pyrat.atm.nlayers, temp, pressure/pc.bar), indent=2)

      # Calculate extinction-coefficient in C:
      extinct_coeff[:] = 0.0
      ec.extinction(extinct_coeff,
          pyrat.voigt.profile, pyrat.voigt.size, pyrat.voigt.index,
          pyrat.voigt.lorentz, pyrat.voigt.doppler,
          pyrat.spec.wn, pyrat.spec.own, pyrat.spec.odivisors,
          density, molq, pyrat.mol.radius, pyrat.mol.mass,
          pyrat.iso.imol, pyrat.iso.mass, pyrat.iso.ratio,
          ziso, pyrat.iso.iext,
          pyrat.lt.wn, pyrat.lt.elow, pyrat.lt.gf, pyrat.lt.isoid,
          pyrat.ex.ethresh, pressure, temp,
          verb-10, int(add),
          int(pyrat.spec.resolution is not None))
      # Store output:
      if grid:   # Into grid
          pyrat.ex.etable[:, itemp, ilayer] = extinct_coeff
      elif add:  # Into ex.ec array for atmosphere
          pyrat.ex.ec[ilayer:ilayer+1] = extinct_coeff
      else:      # return single-layer EC of given layer
          return extinct_coeff


def get_ec(pyrat, layer):
  """
  Compute per-species extinction coefficient at requested layer.
  """
  # Interpolate:
  if pyrat.ex.extfile is not None:
      exc    = np.zeros((pyrat.ex.nmol, pyrat.spec.nwave))
      label  = []
      temp   = pyrat.atm.temp[layer]
      itemp  = np.where(pyrat.ex.temp <= temp)[0][-1]
      if itemp == len(pyrat.ex.temp):
          itemp -= 1
      for i in np.arange(pyrat.ex.nmol):
          imol = np.where(pyrat.mol.ID == pyrat.ex.molID[i])[0][0]
          label.append(pyrat.mol.name[imol])
          etable = pyrat.ex.etable[i,:,layer,:]
          exc[i] = ((etable[itemp  ] * (pyrat.ex.temp[itemp+1] - temp) +
                     etable[itemp+1] * (temp - pyrat.ex.temp[itemp]  ) ) /
                    (pyrat.ex.temp[itemp+1] - pyrat.ex.temp[itemp])      )
          exc[i] *= pyrat.atm.d[layer, imol]
  # Line-by-line:
  else:
      exc = extinction(pyrat, [layer], grid=False, add=False)
      label = []
      for i in np.arange(pyrat.ex.nmol):
          imol = np.where(pyrat.mol.ID == pyrat.ex.molID[i])[0][0]
          exc[i] *= pyrat.atm.d[layer,imol]
          label.append(pyrat.mol.name[imol])
  return exc, label
