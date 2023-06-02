# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Extinction',
    'Optdepth',
    'Cloud',
    'Physics',
]

import numpy as np

from .. import tools as pt
from .. import constants as pc
from .. import atmosphere as pa


class Extinction():
    def __init__(self, inputs, log):
        self.ec = None # line-transition extinction coefficient in cm-1
        self.extfile = inputs.extfile
        # Temperature sampling to compute opacity grid
        self.tmin = inputs.tmin
        self.tmax = inputs.tmax
        self.tstep = inputs.tstep


    def __str__(self):
        fmt = {'float': '{:.2e}'.format}
        fw = pt.Formatted_Write()
        fw.write('Extinction-coefficient information:')
        if self.ec is not None:
            fw.write(
                '\nLBL extinction coefficient for the atmospheric model '
                '(ec, cm-1) [layer, wave]:\n{}', self.ec, fmt=fmt)
        extfile = ['None'] if self.extfile is None else self.extfile
        fw.write("Extinction-coefficient table filename(s) (extfile): {}",
            '\n    '.join(extfile))
        if self.extfile is None:
            return fw.text
        fw.write('Minimum temperature (tmin, K): {:6.1f}', self.tmin)
        fw.write('Maximum temperature (tmax, K): {:6.1f}', self.tmax)
        fw.write('Temperature sampling interval (tstep, K): {self.tstep:6.1f}')
        return fw.text


class Cloud(object):
    """Interface to collect all Cloud opacity models"""
    def __init__(self, model_names, pars, fpatchy=None, log=None):
        self.model_names = model_names
        self.models = []    # List of cloud models
        self.pnames = []
        self.texnames = []
        self.npars = 0
        self.pars = None
        self.ec = None  # Cloud extinction coefficient
        self.fpatchy = fpatchy

        if model_names is None:
            return
        for name in model_names:
            model = pa.clouds.get_model(name)
            self.models.append(model)
            self.npars += model.npars
            self.pnames += model.pnames
            self.texnames += model.texnames

        # Parse parameters:
        if pars is None:
            return
        self.pars = pars
        input_npars = len(self.pars)
        if self.npars != input_npars:
            log.error(
                f'Number of input cloud parameters ({input_npars}) '
                'does not match the number of required '
                f'model parameters ({self.npars})'
            )
        j = 0
        for model in self.models:
            model.pars = self.pars[j:j+model.npars]
            j += model.npars

    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('Cloud-opacity models (models):')
        for model in self.models:
            fw.write('\n' + str(model))
        fw.write('\nPatchiness fraction (fpatchy): {:.3f}', self.fpatchy)
        fw.write('Total atmospheric cloud extinction-coefficient '
                 '(ec, cm-1):\n{}', self.ec, fmt={'float':' {:.3e}'.format})
        return fw.text


class Optdepth():
  def __init__(self, inputs, log):
      self.ec = None  # Total extinction coefficient [nlayers, nwave]
      self.raypath = []  # Distance along ray path  [nlayers]
      self.depth = None  # Optical depth at raypath [nlayers, nwave]
      self.B = None  # Blackbody Planck emission [nlayers, nwave]
      self.ideep = None  # Layer index where depth reached maxdepth [nwave]

      self.maxdepth = inputs.maxdepth  # Maximum optical depth to calculate
      self.rt_path = inputs.rt_path  # Radiative-transfer observing geometry

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Optical depth information:')
      fw.write('Observing geometry (rt_path): {}', self.rt_path)
      if self.ec is not None:
          fw.write(
              'Total atmospheric extinction coefficient (ec, cm-1) [layer'
              ', wave]:\n{}',
              self.ec,
              fmt={'float':'{: .3e}'.format},
          )
      if self.depth is None:
          fw.write(
              '\nMaximum optical depth to calculate (maxdepth): {:.2f}',
              self.maxdepth,
          )
          return fw.text

      ideepest = np.amax(self.ideep)
      if self.rt_path in pc.transmission_rt:
          fw.write(
              '\nDistance along the ray path across each layer '
              '(outside-in) at each impact parameter (raypath, km):',
          )
          with np.printoptions(formatter={'float':'{:.1f}'.format},threshold=6):
              fw.write('    IP[{:3d}]: {}', 1, self.raypath[1]/pc.km)
              fw.write('    IP[{:3d}]: {}', 2, self.raypath[2]/pc.km)
              fw.write('    IP[{:3d}]: {}', 3, self.raypath[3]/pc.km)
              fw.write('    ...')
              fw.write('    IP[{:3d}]: {}', len(self.raypath)-1,
                       self.raypath[len(self.raypath)-1]/pc.km)
          od_text = (
              '\nOptical depth at each impact parameter, down to '
              'max(ideep) (depth):'
          )
      elif self.rt_path in pc.emission_rt:
          fw.write(
              '\nDistance across each layer along a normal ray path '
              '(raypath, km):\n    {}',
              self.raypath/pc.km,
              fmt={'float':'{:.1f}'.format},
              edge=4,
          )
          od_text = (
              '\nOptical depth at each layer along a normal ray '
              'path into the planet, down to max(ideep) (depth):'
          )

      fw.write(
          '\nMaximum optical depth to calculate (maxdepth): {:.2f}',
          self.maxdepth,
      )
      fw.write(
          'Layer index where the optical depth reaches maxdepth (ideep):'
          '\n    {}',
          self.ideep,
          fmt={'int': '{:3d}'.format},
          edge=7,
      )
      fw.write('Maximum ideep (deepest layer reaching maxdepth): {}', ideepest)

      if self.rt_path in pc.emission_rt:
          fw.write(
              '\nPlanck emission down to max(ideep) (B, erg s-1 cm-2 '
              'sr-1 cm):\n{}',
              self.B[0:ideepest+1],
              fmt={'float':'{: .3e}'.format},
          )

      fw.write(
          '{}\n{}', od_text, self.depth[0:ideepest+1],
          fmt={'float':'{: .3e}'.format},
      )
      return fw.text


class Physics(object):
    """Physical properties about the planet and star"""
    def __init__(self, inputs):
        # Stellar properties
        self.tstar = inputs.tstar
        self.rstar = inputs.rstar
        self.mstar = inputs.mstar
        self.log_gstar = inputs.log_gstar

        # Stellar spectrum filename
        self.starspec = inputs.starspec
        self.kurucz = inputs.kurucz
        self.marcs = inputs.marcs
        self.phoenix = inputs.phoenix

        self.starwn = None  # Input stellar wavenumber array
        self.starflux = None  # Input stellar flux spectrum in  FINDME units

    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('Physical properties information:')

        rstar = pt.none_div(self.rstar, pc.rsun)
        mstar = pt.none_div(self.mstar, pc.msun)
        fw.write(
            '\nStellar effective temperature (tstar, K): {:.1f}',
            self.tstar,
        )
        fw.write('Stellar radius (rstar, Rsun): {:.3f}', rstar)
        fw.write('Stellar mass (mstar, Msun):   {:.3f}', mstar)
        fw.write(
            'Stellar surface gravity (log_gstar, cm s-2): {:.2f}',
            self.log_gstar,
        )
        #fw.write('Planet-to-star radius ratio (rprs):   {:.5f}', rprs)
        if self.starspec is not None:
            fw.write(f"Input stellar spectrum (starspec): '{self.starspec}'")
        elif self.kurucz is not None:
            fw.write(f"Input Kurucz stellar spectrum (kurucz): '{self.kurucz}'")
        elif self.marcs is not None:
            fw.write(f"Input MARCS stellar spectrum (marcs): '{self.marcs}'")
        elif self.phoenix is not None:
            fw.write(
                f"Input PHOENIX stellar spectrum (phoenix): '{self.phoenix}'",
            )
        elif self.starflux is not None:
            fw.write(
                "Input stellar spectrum is a blackbody at Teff = {:.1f} K.",
                self.tstar,
            )
        fw.write('Stellar spectrum wavenumber (starwn, cm-1):\n    {}',
            self.starwn, fmt={'float': '{:10.3f}'.format})
        fw.write('Stellar flux spectrum (starflux, erg s-1 cm-2 cm):\n    {}',
            self.starflux, fmt={'float': '{: .3e}'.format})
        return fw.text
