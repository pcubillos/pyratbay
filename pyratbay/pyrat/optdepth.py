import sys, os
import time
import numpy as np
import scipy.integrate as si

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

  pt.msg(pyrat.verb, "\nBegin optical-depth calculation.")
  ti = time.time()
  # Flag to indicate that the extinction has been computed at given layer:
  computed = np.zeros(pyrat.atm.nlayers, np.short)

  # Evaluate the extinction coefficient at each layer:
  pyrat.ex.ec    = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
  pyrat.od.ec    = np.empty((pyrat.atm.nlayers, pyrat.spec.nwave))
  pyrat.od.depth = np.empty((pyrat.atm.nlayers, pyrat.spec.nwave))
  pyrat.od.ideep = np.tile(pyrat.atm.nlayers-1, pyrat.spec.nwave)
  print("Init:   {:.6f}".format(time.time()-ti))

  # Calculate the ray path:
  ti = time.time()
  path(pyrat, pyrat.atm.radius, pyrat.od.path)
  print("Path:   {:.6f}".format(time.time()-ti))

  # Obtain the extinction-coefficient:
  ti = time.time()
  r = 0
  while r < pyrat.atm.nlayers:
    # Interpolate from table:
    if pyrat.ex.extfile is not None:
      ec.interp_ec(pyrat.ex.ec[r],
                   pyrat.ex.etable[:,:,r,:], pyrat.ex.temp, pyrat.ex.molID,
                   pyrat.atm.temp[r], pyrat.atm.d[r], pyrat.mol.ID)
    # On-the-spot calculation of the extinction coefficient:
    elif pyrat.lt.nTLI > 0:
      ex.extinction(pyrat, pyrat.ex.ec[r:r+1], r,
                    pyrat.atm.temp[r], pyrat.iso.z[:,r], add=1)
    r += 1
  print("Interp: {:.6f}".format(time.time()-ti))


  ti = time.time()
  r = 0
  while r < pyrat.atm.nlayers:
    # Sum all contributions to the extinction:
    pyrat.od.ec[r] = pyrat.ex.ec[r] + pyrat.cs.ec[r] + pyrat.haze.ec[r]
    r += 1
  print("Add:    {:.6f}".format(time.time()-ti))


  ti = time.time()
  # Calculate the optical depth for each wavenumber:
  i = 0
  if pyrat.od.path == "eclipse":
    while i < pyrat.spec.nwave:
      pyrat.od.ideep[i] = t.cumtrapz(pyrat.od.depth[:,i], pyrat.od.ec[:,i],
                                     pyrat.od.raypath, pyrat.od.maxdepth)
      i += 1
  else: # pyrat.od.path == "transit"
    while i < pyrat.spec.nwave:
      r = 0
      while r < pyrat.atm.nlayers:
        # Optical depth at each level (tau = integral e*ds):
        pyrat.od.depth[r,i] = t.trapz(pyrat.od.ec[:r+1,i], pyrat.od.raypath[r])

        # Stop calculating the op. depth at this wavenumber if reached maxdeph:
        if pyrat.od.depth[r,i] >= pyrat.od.maxdepth:
          pyrat.od.ideep[i] = r
          break
        r += 1
      i += 1
  print("Integ:  {:.6f}".format(time.time()-ti))
  pt.msg(pyrat.verb, "Done.")


def path(pyrat, radius, path):
  """
  Calculate the distance along the ray path over each interval (layer).

  Notes:
  ------
  - Note that for eclipse geometry the path is always the same.  However,
    for transit geometry the path is unique to each impact parameter; hence,
    the calculation is more laborious.

  Parameters:
  -----------
  pyrat: Pyrat instance
  radius: 1D float ndarray
  path: String
     Flag that indicates eclipse or transit geometry.
  """
  if   path == "eclipse":
    diffrad = np.empty(pyrat.atm.nlayers-1, np.double)
    cu.ediff(radius, diffrad, pyrat.atm.nlayers)
    pyrat.od.raypath = -diffrad

  elif path == "transit":
    pyrat.od.raypath = []
    r = 0
    # Compute the path for each impact parameter:
    while r < pyrat.atm.nlayers:
      raypath = np.empty(r, np.double)
      for i in np.arange(r):
        raypath[i] = (np.sqrt(radius[i  ]**2 - radius[r]**2) -
                      np.sqrt(radius[i+1]**2 - radius[r]**2) )
      pyrat.od.raypath.append(raypath)
      pt.msg(pyrat.verb-10,
        "Raypath[{:3d}]: {}".format(r, pyrat.od.raypath[r]), 2)
      r += 1

  return

