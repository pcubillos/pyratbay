# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["specpercent"]

import sys
import os
if sys.version_info.major == 3:
    import configparser
else:
    import ConfigParser as configparser
import ctypes
import multiprocessing as mp

import numpy as np

from .. import constants as pc
from .. import starspec  as ps
from .. import io        as io

sys.path.append(pc.ROOT + "modules/MCcubed/")
import MCcubed as mc3


def worker(pyrat, idx, models, post, uind):
  """
  Multiprocessing worker that computes the model spectra.

  Parameters
  ----------
  pyrat: A Pyrat object
  idx: List of integers
     Indices in uind to draw (also the indices in models).
  models: 2D float ndarray
     Shared memory array where to store the models (of shape [nunique,nwave]).
  post: 2D float ndarray
     Posterior sample of parameters to compute (if shape [nsamples,npars]).
  uind: 1D integer ndarray
     Indices of unique samples in post (of shape nunique).
  """
  params = pyrat.ret.params
  ifree = np.where(pyrat.ret.stepsize >0)[0]
  j = 0
  for i in idx:
    params[ifree] = post[uind[i]]
    models[i], dummy = pyrat.eval(params)
    j += 1
    if 0 in idx and j%(len(idx)/10)==0:
      print("{:3.0f}% completed.".format(j*100.0/len(idx)))


def specpercent(cfile, ncpu=None, nmax=None):
  """
  Compute the upper and lower 1sigma and 2sigma regions of an
  MCMC sampled spectra.  Sigma boundaries are calculated from the
  0.683 and 0.954 percentiles of the distribution at each wavelength.

  Parameters
  ----------
  cfile: String
     A Pyrat-bay MCMC configuration file.
  ncpu: Integer
     Number of parallel CPUs to use (defaults to the maximum available).
  nmax: Integer
     Maximum number of samples to include (defaults to the size of the
     sample minus the originally burned-in samples).

  Returns
  -------
  wl: 1D float ndarray
     Spectrum wavelength array (in microns).
  bestmodel: 1D float ndarray
     Best-fitting models of the posterior sample.
  low2: 1D float ndarray
     Spectrum values of the lower 0.954 percentile of the posterior.
  low1: 1D float ndarray
     Spectrum values of the lower 0.683 percentile of the posterior.
  high1: 1D float ndarray
     Spectrum values of the upper 0.683 percentile of the posterior.
  high2: 1D float ndarray
     Spectrum values of the upper 0.954 percentile of the posterior.
  """
  # Need to put the import here to avoid circular imports:
  from .. import Pyrat

  # Set number of cpus:
  if ncpu is None:
      ncpu = mp.cpu_count()  # Default to max
  else:
      ncpu = np.clip(ncpu, 1, mp.cpu_count())

  # Initialize Pyrat object from given config file:
  path, cfg = os.path.split(os.path.realpath(cfile))
  path += "/"

  log = mc3.utils.Log(None, width=80)  # Avoid overwriting the log file
  pyrat = Pyrat(cfile, log=log)
  wl    = 1.0/(pyrat.spec.wn*pc.um)
  nwave = len(wl)
  pyrat.verb = -5

  # Extract MCMC setup info:
  config = configparser.ConfigParser()
  config.read([cfile])
  rootname = config.get('pyrat','logfile').replace('.log','')
  savefile = rootname + '.npz'
  burn = (config.getint('pyrat','burnin')
         *config.getint('pyrat','nchains')
         /config.getint('pyrat','thinning'))

  # Burned-in posterior sample:
  post = np.load(path + savefile)["Z"]
  nsamples = np.shape(post)[0]
  if nmax is None:
    nmax = nsamples  # Default to the size of the sample
  ini  = nsamples-nmax
  # Take maximum between burn index and user-input index:
  burn = int(np.amax([burn,ini]))
  post = post[burn:]

  # Unique MCMC samples:
  u, uind, uinv = np.unique(post[:,0], return_index=True, return_inverse=True)
  nruns = len(u)
  print("\nComputing {:d} models:".format(nruns))

  # Shared memory array:
  sm_models = mp.Array(ctypes.c_double, nruns*nwave)
  models = np.ctypeslib.as_array(sm_models.get_obj()).reshape((nruns,nwave))

  # Compute spectra:
  processes = []
  for n in np.arange(ncpu):
    start =  n    * nruns//ncpu
    end   = (n+1) * nruns//ncpu
    idx = np.arange(start, end)
    proc = mp.Process(target=worker, args=(pyrat, idx, models, post, uind))
    processes.append(proc)
    proc.start()
  # Make sure all processes finish their work:
  for n in np.arange(ncpu):
    processes[n].join()

  print("\nComputing percentiles:")
  # Compute 1sigma, 2sigma percentiles:
  #   1 sigma: 0.6827 -- 0.15865
  #   2 sigma: 0.9545 -- 0.02275
  low1    = np.zeros(nwave)
  low2    = np.zeros(nwave)
  high1   = np.zeros(nwave)
  high2   = np.zeros(nwave)
  for i in np.arange(nwave):
    msample = models[uinv,i]
    low2[i]  = np.percentile(msample,  2.275)
    low1[i]  = np.percentile(msample, 15.865)
    high1[i] = np.percentile(msample, 100-15.865)
    high2[i] = np.percentile(msample, 100- 2.275)
    if i % int(nwave//10) == 0:
      print("{:3.0f}% completed.".format(i*100.0/nwave))

  # Load best-fitting model:
  wn, bestmodel = io.read_pyrat(path + rootname + "_bestfit_spectrum.dat")

  # Write to file:
  with open(path + rootname + "_percentiles.dat", "w") as f:
    f.write("# wl (um)    Depth        -2sigma      -1sigma      "
                         "+1sigma      +2sigma\n")
    for i in np.arange(nwave):
      f.write("{:.5e}  {:.5e}  {:.5e}  {:.5e}  {:.5e}  {:.5e}\n".
          format(wl[i], bestmodel[i], low2[i], low1[i], high1[i], high2[i]))

  log.close()
  os.remove("dummy.log")
  return wl, bestmodel, low2, low1, high1, high2
