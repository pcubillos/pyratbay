# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys
import ctypes
import multiprocessing as mpr

import numpy as np

from .  import extinction as ex
from .. import constants  as pc

sys.path.append(pc.ROOT + 'lib')
import extcoeff as ec
import cutils   as cu
import trapz    as t


def opticaldepth(pyrat):
  """
  Calculate the optical depth.
  """
  od = pyrat.od
  pyrat.log.msg('\nBegin optical-depth calculation.')

  # Evaluate the extinction coefficient at each layer:
  pyrat.ex.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
  od.ec       = np.empty((pyrat.atm.nlayers, pyrat.spec.nwave))
  od.depth    = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
  od.ideep    = np.tile(pyrat.atm.nlayers-1, pyrat.spec.nwave)
  if pyrat.cloud.fpatchy:
      od.epatchy = np.empty((pyrat.atm.nlayers, pyrat.spec.nwave))
      od.pdepth  = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  rtop = pyrat.atm.rtop
  # Calculate the ray path:
  path(pyrat)

  # Interpolate extinction coefficient from table:
  if pyrat.ex.extfile is not None:
      r = rtop
      while r < pyrat.atm.nlayers:
          ec.interp_ec(pyrat.ex.ec[r],
                       pyrat.ex.etable[:,:,r,:], pyrat.ex.temp, pyrat.ex.molID,
                       pyrat.atm.temp[r], pyrat.atm.d[r], pyrat.mol.ID)
          r += 1

  # Calculate the extinction coefficient on the spot:
  elif pyrat.lt.tlifile is not None:
      sm_ext = mpr.Array(ctypes.c_double,
          np.zeros(pyrat.atm.nlayers*pyrat.spec.nwave, np.double))
      pyrat.ex.ec = np.ctypeslib.as_array(sm_ext.get_obj()).reshape(
          (pyrat.atm.nlayers, pyrat.spec.nwave))
      processes = []
      indices = np.arange(rtop, pyrat.atm.nlayers) % pyrat.ncpu
      for i in range(pyrat.ncpu):
          proc = mpr.Process(target=ex.extinction,   #      grid   add
                      args=(pyrat, np.where(indices==i)[0], False, True))
          processes.append(proc)
          proc.start()
      for proc in processes:
          proc.join()

  # Sum all contributions to the extinction (except clouds):
  od.ec[rtop:] = (pyrat.ex.ec      [rtop:] +
                  pyrat.cs.ec      [rtop:] +
                  pyrat.rayleigh.ec[rtop:] +
                  pyrat.alkali.ec  [rtop:])
  # Add cloud if not fpatchy, else separate Eclear and Ecloudy:
  if pyrat.cloud.fpatchy is None:
      od.ec[rtop:] += pyrat.cloud.ec[rtop:]
  else:
      od.epatchy[rtop:] = np.copy(od.ec[rtop:]) + pyrat.cloud.ec[rtop:]

  # Calculate the optical depth for each wavenumber:
  if od.path == 'eclipse':
      i = 0
      while i < pyrat.spec.nwave:
          od.ideep[i] = t.cumtrapz(od.depth  [rtop:,i],
                                   od.ec     [rtop:,i],
                                   od.raypath[rtop:],
                                   od.maxdepth) + rtop
          i += 1
  elif od.path == 'transit':
      od.ideep = np.array(od.ideep, dtype=np.intc)
      od.ideep[:] = -1
      r = rtop
      while r < pyrat.atm.nlayers:
          # Optical depth at each level (tau = integral e*ds):
          od.depth[r] = t.optdepth(od.ec[rtop:r+1], od.raypath[r],
                                   od.maxdepth, od.ideep, r)
          if pyrat.cloud.fpatchy is not None:
              od.pdepth[r] = t.optdepth(od.epatchy[rtop:r+1], od.raypath[r],
                                        np.inf, od.ideep, r)
          r += 1
      od.ideep[od.ideep<0] = pyrat.atm.nlayers - 1
  pyrat.log.msg('Optical depth done.')


def path(pyrat):
  """
  Calculate the distance along the ray path over each interval (layer).

  Notes
  -----
  For eclipse geometry, the path is always the same.
  For transit geometry, the path is unique to each impact parameter.
  """
  if pyrat.od.path == 'eclipse':
      radius = pyrat.atm.radius
      diffrad = np.empty(pyrat.atm.nlayers-1, np.double)
      cu.ediff(radius, diffrad, pyrat.atm.nlayers)
      pyrat.od.raypath = -diffrad

  elif pyrat.od.path == 'transit':
      pyrat.od.raypath = []
      radius  = pyrat.atm.radius[pyrat.atm.rtop:]
      nlayers = pyrat.atm.nlayers - pyrat.atm.rtop
      # Empty-filling layers that don't contribute:
      for r in range(pyrat.atm.rtop):
          pyrat.od.raypath.append([])
      # Compute the path for each impact parameter:
      r = 0
      while r < nlayers:
          raypath = np.empty(r, np.double)
          for i in range(r):
              raypath[i] = (np.sqrt(radius[i  ]**2 - radius[r]**2) -
                            np.sqrt(radius[i+1]**2 - radius[r]**2) )
          pyrat.od.raypath.append(raypath)
          pyrat.log.msg('Raypath[{:3d}]: {}'.format(r, pyrat.od.raypath[r]),
                        verb=4, indent=2)
          r += 1

