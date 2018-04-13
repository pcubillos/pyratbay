# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import numpy as np
import scipy.special as ss
from scipy.ndimage.filters import gaussian_filter1d as gaussf

from .. import pyrat      as py
from .. import atmosphere as pa
from .. import constants  as pc
from .. import blackbody  as bb


rootdir = os.path.realpath(os.path.dirname(__file__) + "/../../")
TEAdir = rootdir + "/modules/TEA/"


def radeq(args):
  """
  Compute radiative-thermochemical equilibrium atmosphere.

  Parameters
  ----------
  args: Dict
     PyratBay arguments dictionary

  Returns
  -------
  pyrat: Pyrat object
     The Pyrat object at the last iteration.
  temp: 2D float ndarray
     Temperature profiles (in Kelvin) at each iteration.
  Fup: 2D float ndarray
     Upwards flux spectrum through each layer under the two-stream
     approximation at the last iteration (erg s-1 cm-2 cm).
  Fdown: 2D float ndarray
     Downwards flux spectrum through each layer under the two-stream
     approximation at the last iteration (erg s-1 cm-2 cm).
  Qup: 2D float ndarray
     Upwards bolometric net flux through each layer under the
     two-stream approximation at each iteration (erg s-1 cm-2).
  Qdown: 2D float ndarray
     Downwards bolometric net flux through each layer under the
     two-stream approximation at each iteration (erg s-1 cm-2).
  """
  # Load Pyrat object:
  pyrat = py.init(args.cfile)

  pyrat.verb = 0  # Mute it
  pyrat.od.maxdepth = np.inf
  wn      = pyrat.spec.wn
  nlayers = pyrat.atm.nlayers
  tmin = np.amax((pyrat.cs.tmin, pyrat.ex.tmin))
  tmax = np.amin((pyrat.cs.tmax, pyrat.ex.tmax))
  elements, stoich = pa.stoich(pyrat.mol.name)

  # Default inputs:
  TEAdc   = 5      # TEA duty cycle
  fscale0 = 1.0e5  # Initial temperature scale factor
  maxf    = 2.0e5  # Maximum temperature scale factor

  # Top boundary condition:
  Fstar = args.beta * (0.5*pyrat.phy.rstar/pyrat.phy.smaxis)**2 \
          * pyrat.phy.starflux

  # Set internal net bolometric flux to sigma*Tint^4:
  Fint = pc.sigma*pyrat.phy.tint**4 \
         * bb.Bwn(wn, pyrat.phy.tint) / np.trapz(bb.Bwn(wn,pyrat.phy.tint),wn)

  # Arrays to iterate over:
  temp  = np.zeros((args.nsamples+1, nlayers))
  Qup   = np.zeros((args.nsamples+1, nlayers))
  Qdown = np.zeros((args.nsamples+1, nlayers))
  dF    = np.zeros((args.nsamples+1, nlayers))
  fscale = np.tile(fscale0, (args.nsamples+1,nlayers))

  Fsign = np.tile(2.0, nlayers)
  Fdown = np.zeros((nlayers, pyrat.spec.nwave))
  Fup   = np.zeros((nlayers, pyrat.spec.nwave))
  sigma = 1.0

  # Starting point:
  temp[0]  = pyrat.atm.temp
  Fdown[0] = Fstar
  q = np.copy(pyrat.atm.q)

  print("Radiative-thermochemical equilibrium calculation:")
  for k in np.arange(args.nsamples):
    sys.stdout.write("\rIteration {:3d}/{:d}.".format(k+1, args.nsamples))
    sys.stdout.flush()
    # Update abundances every TEAdc iterations:
    if ((k+1) % TEAdc) == 0:
      q = pa.calcatm(args.atmfile, pyrat.atm.press, temp[k], pyrat.mol.name,
                     xsolar=args.xsolar, elements=elements, nproc=args.nproc)
    # Update opacities:
    pyrat = py.run(pyrat, [temp[k], q])
    B = pyrat.od.B

    # Diffusivity factor (Eqs 147 and B5 of Heng2014):
    dtau0 = np.diff(pyrat.od.depth, n=1, axis=0)
    D = -1/dtau0 * np.log((1-dtau0)*np.exp(-dtau0) + dtau0**2*ss.expn(1,dtau0))
    D[~np.isfinite(D)] = 1.0
    trans = np.exp(-D*dtau0)
    Bp = np.diff(B, n=1, axis=0) / dtau0

    # Diffuse approximation to compute downward and upward fluxes:
    for i in np.arange(nlayers-1):
      Fdown[i+1] = trans[i]*Fdown[i] \
               + np.pi*(B[i+1] - trans[i]*B[i] - 0.5*Bp[i]*(1-trans[i]))
    Fup[nlayers-1] = Fdown[nlayers-1] + Fint
    for i in np.arange(nlayers-1)[::-1]:
      Fup[i] = trans[i]*Fup[i+1] \
               + np.pi*(B[i] - trans[i]*B[i+1] + 0.5*Bp[i]*(1-trans[i]))

    # Bolometric fluxes:
    for i in np.arange(nlayers):
      Qup  [k,i] = np.trapz(Fup[i],   wn)
      Qdown[k,i] = np.trapz(Fdown[i], wn)

    # Bolometric net flux through each layer:
    dQ = Qup[k] - Qdown[k]
    dF[k] = np.ediff1d(dQ, to_begin=0)

    # Update scaling factor:
    idiff = np.sign(dF[k]) != Fsign
    fscale[k+1,  idiff] = fscale[k,  idiff]*0.9
    fscale[k+1, ~idiff] = fscale[k, ~idiff]*1.05
    fscale[k+1] = np.clip(fscale[k+1], 0, maxf)
    Fsign = np.sign(dF[k])

    # Adaptive temperature update:
    dT = fscale[k+1]*Fsign*np.abs(dF[k])**0.1 / \
          (np.ediff1d(np.log(pyrat.atm.press), to_begin=1.0)
           * pc.sigma*pyrat.atm.temp**3)
    temp[k+1] = pyrat.atm.temp + dT
    temp[k+1,0] = temp[k+1,1]  # Isothermal top
    temp[k+1] = np.clip(temp[k+1], tmin, tmax)

    # Smooth out kinks:
    dtm = np.mean(np.abs(np.diff(temp[:k+1],n=1,axis=0)),axis=1)
    if np.size(dtm) > 0:
      sigma = np.clip(0.5*np.mean(dtm[-5:]), 0.5, 1.0)
    temp[k+1] = gaussf(temp[k+1], sigma)

  # Update last temp iteration:
  pyrat.atm.temp = temp[k+1]

  return pyrat, temp, Fup, Fdown, Qup, Qdown
