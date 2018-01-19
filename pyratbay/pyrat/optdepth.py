# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import time
import ctypes
import numpy as np
import scipy.integrate as si
import multiprocessing as mpr

from .. import tools as pt

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../lib')
import extinction as ex
import extcoeff   as ec
import cutils     as cu
import trapz      as t


def opticaldepth(pyrat):
  """
  Calculate the optical depth.
  """

  pt.msg(pyrat.verb-3, "\nBegin optical-depth calculation.", pyrat.log)
  ti = time.time()
  # Flag to indicate that the extinction has been computed at given layer:
  computed = np.zeros(pyrat.atm.nlayers, np.short)

  # Evaluate the extinction coefficient at each layer:
  pyrat.ex.ec    = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
  pyrat.od.ec    = np.empty((pyrat.atm.nlayers, pyrat.spec.nwave))
  pyrat.od.depth = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
  pyrat.od.ideep = np.tile(pyrat.atm.nlayers-1, pyrat.spec.nwave)
  if pyrat.haze.fpatchy:
    pyrat.od.epatchy = np.empty((pyrat.atm.nlayers, pyrat.spec.nwave))
    pyrat.od.pdepth  = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
  #print("Init:   {:.6f}".format(time.time()-ti))

  rtop = pyrat.atm.rtop
  # Calculate the ray path:
  ti = time.time()
  path(pyrat)
  #print("Path:   {:.6f}".format(time.time()-ti))

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
    indices = np.arange(rtop, pyrat.atm.nlayers) % pyrat.nproc
    for i in np.arange(pyrat.nproc):
      proc = mpr.Process(target=ex.extinction,   #      grid   add
                  args=(pyrat, np.where(indices==i)[0], False, True))
      processes.append(proc)
      proc.start()
    for i in np.arange(pyrat.nproc):
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

  ti = time.time()
  # Calculate the optical depth for each wavenumber:
  i = 0
  if pyrat.od.path == "eclipse":
    while i < pyrat.spec.nwave:
      pyrat.od.ideep[i] = t.cumtrapz(pyrat.od.depth  [rtop:,i],
                                     pyrat.od.ec     [rtop:,i],
                                     pyrat.od.raypath[rtop:],
                                     pyrat.od.maxdepth) + rtop
      i += 1
  else: # pyrat.od.path == "transit"
    while i < pyrat.spec.nwave:
      r = rtop
      while r < pyrat.atm.nlayers:
        # Optical depth at each level (tau = integral e*ds):
        pyrat.od.depth[r,i] = t.trapz(pyrat.od.ec[rtop:r+1,i],
                                      pyrat.od.raypath[r])
        if pyrat.haze.fpatchy is not None:
          pyrat.od.pdepth[r,i] = t.trapz(pyrat.od.epatchy[rtop:r+1,i],
                                         pyrat.od.raypath[r])

        # Stop calculating the op. depth at this wavenumber if reached maxdeph:
        if pyrat.od.depth[r,i] >= pyrat.od.maxdepth:
          pyrat.od.ideep[i] = r
          break
        r += 1
      i += 1
  #print("Integ:  {:.6f}".format(time.time()-ti))
  pt.msg(pyrat.verb-3, "Done.", pyrat.log)


def path(pyrat):
  """
  Calculate the distance along the ray path over each interval (layer).

  Notes:
  ------
  - Note that for eclipse geometry the path is always the same.  However,
    for transit geometry the path is unique to each impact parameter; hence,
    the calculation is more laborious.
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
      pt.msg(pyrat.verb-6, "Raypath[{:3d}]: {}".
                            format(r, pyrat.od.raypath[r]), pyrat.log, 2)
      r += 1

  return

