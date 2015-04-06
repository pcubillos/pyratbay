import numpy as np
import scipy.integrate as si
import time

import extinction as ex
import extcoeff   as ec
import ptools     as pt
import simpson    as s

def opticaldepth(pyrat):
  """
  Calculate the optical depth.
  """

  # Allocate tau, last, intens:

  # Flag to indicate that the extinction has been computed at given layer:
  computed = np.zeros(pyrat.atm.nlayers, np.short)

  # Evaluate the extinction coefficient at each layer:
  pyrat.ex.ec    = np.zeros((pyrat.atm.nlayers, pyrat.nspec))
  pyrat.od.ec    = np.zeros((pyrat.atm.nlayers, pyrat.nspec))
  pyrat.od.depth = np.zeros((pyrat.atm.nlayers, pyrat.nspec))
  pyrat.od.ideep = np.zeros((pyrat.nspec), np.double)

  hsum, hratio, hfactor = [], [], []
  # Calculate the optical depth for each wavenumber:
  for i in np.arange(pyrat.nspec):
    t0 = time.time()
    for r in np.arange(pyrat.atm.nlayers):
      if not computed[r]:
        # Calculate the raypath for this layer:
        path(pyrat, pyrat.atm.radius, r, pyrat.path)
        pt.msg(pyrat.verb-10,
               "Raypath[{:3d}]: {}".format(r, pyrat.od.raypath[r]), 2)
        # Calculate values for Simpson integration:
        hs, hr, hf = s.geth(pyrat.od.raypath[r])
        hsum.append(hs)
        hratio.append(hr)
        hfactor.append(hf)

        # On-the-spot calculation of the extinction coefficient:
        if pyrat.ex.extfile is None:
          #print("FLAG 1, EXTFILE is NONE")
          ex.extinction(pyrat, pyrat.ex.ec[r:r+1], r,
                        pyrat.atm.temp[r], pyrat.iso.z[:,r], add=1)
        # Interpolate from table:
        else:
          #print("FLAG 2, EXTFILE EXISTS")
          ti = time.time()
          ec.interp_ec(pyrat.ex.ec[r],
                       pyrat.ex.etable[:,:,r,:], pyrat.ex.temp, pyrat.ex.molID,
                       pyrat.atm.temp[r], pyrat.atm.d[r], pyrat.mol.ID)
          tf = time.time()
          print("Interpol time: {}\n".format(tf-ti))

        # Sum all contributions to the extinction:
        pyrat.od.ec[r] = pyrat.ex.ec[r] + pyrat.cia.ec[r]
        pt.msg(pyrat.verb-10, "Total EC: {}".format(pyrat.od.ec[:r+1,i]), 2)

        # Update the computed flag:
        computed[r] = 1

      # Optical depth at each level (tau = inttegral e*ds):
      pyrat.od.depth[r,i] = s.simps(pyrat.od.ec[:r+1,i],
                    pyrat.od.raypath[r], hsum[r], hratio[r], hfactor[r])

      if i == 1:
        pt.msg(pyrat.verb-5, "Optical depth at {:7.2f} cm-1 ({:3d}): {:.4e}".
                  format(pyrat.wn[i], r, pyrat.od.depth[r,i]))
        pt.msg(pyrat.verb-5, "Extinction: CIA={:.5e},  Ext={:.5e}".format(
                              pyrat.cia.ec[r,i], pyrat.ex.ec[r,i]))
      pt.msg(pyrat.verb-10, "Optical depth[{:3d}, {:4d}]: {:.4e}".format(
                             r, i, pyrat.od.depth[r,i]), 2)

      # Stop calculating the op. depth at this wavenumber if reached maxdeph:
      if pyrat.od.depth[r,i] >= pyrat.maxdepth:
        pyrat.od.ideep[i] = r
        break

      #print("")
      #pyrat.od.depth[1:,i] = si.cumtrapz(pyrat.ex.ec[:,i], raypath)
    t1 = time.time()
    #print("Optical-depth time: {}\n".format(t1-t0))


def path2(pyrat, radius, r, path):
  """ 
  Calculate the distance along the ray path.

  Parameters:
  -----------
  pyrat: Object
  radius: 1D float ndarray
  r: Integer
  path: String

  """
  if path == "eclipse":
    # Distance traveled perpendiculary inwards:
    raypath = radius[0] - radius

  # Calculate slant path for impact parameters at each layer radius:
  elif path == "transit":
    # Cumulative distance along the raypath from the outermost layer
    # to each layer:
    raypath = (np.sqrt(radius[0    ]**2 - radius[r]**2) -
               np.sqrt(radius[0:r+1]**2 - radius[r]**2) )

  pyrat.od.raypath.append(raypath)
  return


def path(pyrat, radius, r, path):
  """
  Calculate the distance along the ray path over each layer.
  """
  if   path == "eclipse":
    pyrat.od.raypath.append(np.ediff1d(radius))

  elif path == "transit":
    raypath = np.empty(r, np.double)
    for i in np.arange(r):
      raypath[i] = (np.sqrt(radius[i  ]**2 - radius[r]**2) -
                    np.sqrt(radius[i+1]**2 - radius[r]**2) )
    pyrat.od.raypath.append(raypath)

  return

