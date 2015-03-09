import numpy as np
import scipy.integrate as si
import time

import extinction as ex
import extcoeff   as ec

def opticaldepth(pyrat):
  """
  Calculate the optical depth.
  """

  # Allocate tau, last, intens:
  # Define computed:

  # Evaluate the extinction coefficient at each layer:
  pyrat.ex.ec    = np.zeros((pyrat.atm.nlayers, pyrat.nspec))
  pyrat.od.depth = np.zeros((pyrat.atm.nlayers, pyrat.nspec))

  # Calculate the optical depth for each wavenumber:
  for i in np.arange(pyrat.nspec):
    for r in np.arange(pyrat.atm.nlayers):
      if pyrat.ex.extfile is None:
        print("FLAG 1, EXTFILE is NONE")
        # On-the-spot calculation of the extinction coefficient:
        ex.extinction(pyrat, pyrat.ex.ec[r:r+1], r,
                      pyrat.atm.temp[r], pyrat.iso.z[:,r], add=1)
      else:
        print("FLAG 2, EXTFILE EXISTS")
        # Interpolate from table:
        ti = time.time()
        ec.interp_ec(pyrat.ex.ec[r],
                     pyrat.ex.etable[:,:,r,:], pyrat.ex.temp, pyrat.ex.molID,
                     pyrat.atm.temp[r], pyrat.atm.d[r], pyrat.mol.ID)
        tf = time.time()
        print("Interpol time: {}\n".format(tf-ti))

      # Calculate the raypath for this layer:
      raypath = path(pyrat.atm.radius, r, pyrat.path)

      # Optical depth at each level:
      pyrat.od.depth[1:,i] = si.cumtrapz(pyrat.ex.ec[:,i], pyrat.od.raypath)


def path(radius, r, path):
  """ 
  Calculate the distance along the ray path.
  """
  if path == "eclipse":
    # Distance traveled perpendiculary inwards:
    raypath = radius[0] - radius

  # Calculate slant path for impact parameters at each layer radius:
  elif path == "transit":
    # Cumulative distance along the raypath from the outermost layer
    # to each layer:
    raypath = (np.sqrt(radius[0]**2 - radius[r]**2) -
               np.sqrt(radius**2    - radius[r]**2) )

  return raypath
