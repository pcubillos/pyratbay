# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import ctypes
import multiprocessing as mpr

import numpy as np

from .  import extinction as ex
from .. import constants  as pc

sys.path.append(pc.ROOT + 'lib')
import extcoeff   as ec
import cutils     as cu
import trapz      as t


def opticaldepth(pyrat):
  """
  Calculate the optical depth.
  """

  pyrat.log.msg("\nBegin optical-depth calculation.")

  # Evaluate the extinction coefficient at each layer:
  pyrat.ex.ec    = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
  pyrat.od.ec    = np.empty((pyrat.atm.nlayers, pyrat.spec.nwave))
  pyrat.od.depth = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
  pyrat.od.ideep = np.tile(pyrat.atm.nlayers-1, pyrat.spec.nwave)
  if pyrat.haze.fpatchy:
    pyrat.od.epatchy = np.empty((pyrat.atm.nlayers, pyrat.spec.nwave))
    pyrat.od.pdepth  = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  rtop = pyrat.atm.rtop
  # Calculate the ray path:
  path(pyrat)

  # Obtain the extinction-coefficient:
  # Interpolate from table:
  if pyrat.ex.extfile is not None:
    r = rtop
    while r < pyrat.atm.nlayers:
      ec.interp_ec(pyrat.ex.ec[r],
                   pyrat.ex.etable[:,:,r,:], pyrat.ex.temp, pyrat.ex.molID,
                   pyrat.atm.temp[r], pyrat.atm.d[r], pyrat.mol.ID)
      r += 1

  # On-the-spot calculation of the extinction coefficient:
  elif pyrat.lt.nTLI > 0:
    # Put pyrat.ex.ec into shared memory:
    sm_ext = mpr.Array(ctypes.c_double,
                     np.zeros(pyrat.atm.nlayers*pyrat.spec.nwave, np.double))
    pyrat.ex.ec = np.ctypeslib.as_array(sm_ext.get_obj()).reshape(
                                   (pyrat.atm.nlayers, pyrat.spec.nwave))
    # Multi-processing extinction calculation (in C):
    processes = []
    # CPU indices
    indices = np.arange(rtop, pyrat.atm.nlayers) % pyrat.ncpu
    for i in np.arange(pyrat.ncpu):
      proc = mpr.Process(target=ex.extinction,   #      grid   add
                  args=(pyrat, np.where(indices==i)[0], False, True))
      processes.append(proc)
      proc.start()
    for i in np.arange(pyrat.ncpu):
      processes[i].join()

  # Sum all contributions to the extinction (except clouds):
  pyrat.od.ec[rtop:] = (pyrat.ex.ec      [rtop:] +
                        pyrat.cs.ec      [rtop:] +
                        pyrat.rayleigh.ec[rtop:] +
                        pyrat.alkali.ec  [rtop:])
  # Add cloud if not fpatchy, else separate Eclear and Ecloudy:
  if pyrat.haze.fpatchy is None:
    pyrat.od.ec[rtop:] += pyrat.haze.ec[rtop:]
  else:
    pyrat.od.epatchy[rtop:] = np.copy(pyrat.od.ec[rtop:]) + pyrat.haze.ec[rtop:]

  # Calculate the optical depth for each wavenumber:
  if pyrat.od.path == "eclipse":
    i = 0
    while i < pyrat.spec.nwave:
      pyrat.od.ideep[i] = t.cumtrapz(pyrat.od.depth  [rtop:,i],
                                     pyrat.od.ec     [rtop:,i],
                                     pyrat.od.raypath[rtop:],
                                     pyrat.od.maxdepth) + rtop
      i += 1
  else: # pyrat.od.path == "transit"
    pyrat.od.ideep = np.array(pyrat.od.ideep, dtype=np.intc)
    pyrat.od.ideep[:] = -1
    r = rtop
    while r < pyrat.atm.nlayers:
      # Optical depth at each level (tau = integral e*ds):
      pyrat.od.depth[r] = t.optdepth(pyrat.od.ec[rtop:r+1],
                            pyrat.od.raypath[r],
                            pyrat.od.maxdepth, pyrat.od.ideep, r)
      if pyrat.haze.fpatchy is not None:
        pyrat.od.pdepth[r] = t.optdepth(pyrat.od.epatchy[rtop:r+1],
                            pyrat.od.raypath[r],
                            np.inf, pyrat.od.ideep, r)
      r += 1
    pyrat.od.ideep[pyrat.od.ideep<0] = pyrat.atm.nlayers - 1
  pyrat.log.msg("Optical depth done.")


def path(pyrat):
  """
  Calculate the distance along the ray path over each interval (layer).

  Notes
  -----
  For eclipse geometry, the path is always the same.
  For transit geometry, the path is unique to each impact parameter.
  """
  if   pyrat.od.path == "eclipse":
    radius = pyrat.atm.radius
    diffrad = np.empty(pyrat.atm.nlayers-1, np.double)
    cu.ediff(radius, diffrad, pyrat.atm.nlayers)
    pyrat.od.raypath = -diffrad

  elif pyrat.od.path == "transit":
    pyrat.od.raypath = []
    radius  = pyrat.atm.radius[pyrat.atm.rtop:]
    nlayers = pyrat.atm.nlayers - pyrat.atm.rtop
    # Empty-filling layers that don't contribute:
    for r in np.arange(pyrat.atm.rtop):
      pyrat.od.raypath.append([])
    # Compute the path for each impact parameter:
    r = 0
    while r < nlayers:
      raypath = np.empty(r, np.double)
      for i in np.arange(r):
        raypath[i] = (np.sqrt(radius[i  ]**2 - radius[r]**2) -
                      np.sqrt(radius[i+1]**2 - radius[r]**2) )
      pyrat.od.raypath.append(raypath)
      pyrat.log.msg("Raypath[{:3d}]: {}".format(r, pyrat.od.raypath[r]),
                    verb=4, indent=2)
      r += 1

  return
