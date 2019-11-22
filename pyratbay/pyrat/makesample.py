# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np
import scipy.interpolate as sip

from .. import tools     as pt
from .. import constants as pc


def make_wavenumber(pyrat):
  """
  Make the wavenumber sample from user inputs.
  Store all values in CGS units.
  """
  # Alias for Pyrat's Spectrum object:
  spec = pyrat.spec
  log  = pyrat.log

  log.head('\nGenerating wavenumber array.')
  # Low wavenumber boundary:
  if pyrat.inputs.wnlow is None:
      if pyrat.inputs.wlhigh is None:
          log.error('Low wavenumber boundary is undefined.  Either set '
                    'wnlow or wlhigh.')
      else:
          spec.wnlow = 1.0 / spec.wlhigh
  elif pyrat.inputs.wlhigh is not None:
      log.warning('Both wnlow ({:.2e} cm-1) and wlhigh ({:.2e} cm) were '
                  'defined.  Pyrat will take wnlow and ignore wlhigh.'.
                  format(spec.wnlow, spec.wlhigh))

  # High wavenumber boundary:
  if pyrat.inputs.wnhigh is None:
      if pyrat.inputs.wllow is None:
          log.error('High wavenumber boundary is undefined.  Either set '
                    'wnhigh or wllow.')
      else:
          spec.wnhigh = 1.0 / spec.wllow
  elif pyrat.inputs.wllow is not None:
      log.warning('Both wnhigh ({:.2e} cm-1) and wllow ({:.2e} cm) were '
                  'defined.  Pyrat will take wnhigh and ignore wllow.'.
                  format(spec.wnhigh, spec.wllow))

  # Consistency check (wnlow < wnhigh):
  if spec.wnlow > spec.wnhigh:
    log.error('Wavenumber low boundary ({:.1f} cm-1) must be larger '
              'than the high boundary ({:.1f} cm-1).'.
              format(spec.wnlow, spec.wnhigh))

  if spec.wnstep is None:
      log.error('Undefined wavenumber sampling step size (wnstep).')

  if spec.wnosamp is None:
      log.error('Undefined wavenumber oversampling factor (wnosamp).')

  # Set wavelength limits based on the wavenumber limits:
  spec.wlhigh = 1.0 / spec.wnlow
  spec.wllow  = 1.0 / spec.wnhigh

  # Make the wavenumber array:
  spec.wn = np.arange(spec.wnlow, spec.wnhigh, spec.wnstep)

  # Re-set final boundary (stay inside given boundaries):
  if spec.wn[-1] != spec.wnhigh:
      log.warning(
          'Final wavenumber boundary modified from {:.4f} cm-1 (input)\n'
          '                                     to {:.4f} cm-1 (Pyrat).'.
          format(spec.wnhigh, spec.wn[-1]))
  # Set the number of spectral samples:
  spec.nwave = len(spec.wn)

  # Make the fine-sampled (oversampled) wavenumber array:
  spec.ownstep = spec.wnstep / spec.wnosamp
  spec.onwave  = (spec.nwave - 1) *  spec.wnosamp + 1
  spec.own = np.linspace(spec.wn[0], spec.wn[-1], spec.onwave)

  # Get list of divisors:
  spec.odivisors = pt.divisors(spec.wnosamp)

  # User-defined output resolution:
  if spec.resolution is not None:
      f = 0.5 / spec.resolution
      g = (1.0-f) / (1.0+f)
      imax = int(np.ceil(np.log(spec.wnlow/spec.wnhigh) / np.log(g))) + 1
      # These are the wavelength edges of each bin:
      dwn = spec.wnhigh * g**np.arange(imax)
      # Central wavelength of each bin:
      spec.wn = 0.5*(dwn[1:]+dwn[:-1])[::-1]
      # Total number of spectral samples:
      spec.nwave = imax - 1


  # Screen output:
  log.msg(f'Initial wavenumber boundary:  {spec.wnlow:.5e} cm-1  '
          f'({spec.wlhigh/pt.u(spec.wlunits):.3e} {spec.wlunits})', indent=2)
  log.msg(f'Final   wavenumber boundary:  {spec.wnhigh:.5e} cm-1  '
          f'({spec.wllow/pt.u(spec.wlunits):.3e} {spec.wlunits})', indent=2)
  if spec.resolution is None:
      log.msg(f'Wavenumber sampling interval: {spec.wnstep:.2g} cm-1', indent=2)
  else:
      log.msg(f'Spectral resolving power: {spec.resolution:.1f}', indent=2)
  log.msg(f'Wavenumber sample size:      {spec.nwave:8d}\n'
          f'Wavenumber fine-sample size: {spec.onwave:8d}\n', indent=2)
  log.head('Wavenumber sampling done.')


def make_atmprofiles(pyrat):
  """
  Define atmospheric-profile layers sampling.

  Notes
  -----
  There are multiple things going on in this functions, here's a summary:
  - Compute hydrostatic-equilibrium radius profile if necessary.
  - Reset (pressure/radius) boundaries if requested.
  - Resample (pressure/radius/abundance/temperature) profiles if requested.
  - Check whether top of atmosphere crosses the Hill radius.
  - Compute partition-function at layers temperatures.
  """
  log = pyrat.log
  log.head('\nGenerating atmospheric profile sample.')

  # Pyrat and user-input atmospheric-data objects:
  atm    = pyrat.atm
  atm_in = pyrat.inputs.atm

  # Check that the layers are sorted from the top to the bottom of
  #  the atmosphere:
  sort    = np.all(np.ediff1d(atm_in.press) > 0)  # Top to bottom
  reverse = np.all(np.ediff1d(atm_in.press) < 0)  # Bottom to top
  if atm_in.radius is not None:
      sort    *= np.all(np.ediff1d(atm_in.radius) < 0)
      reverse *= np.all(np.ediff1d(atm_in.radius) > 0)

  if sort:       # Layers are in the correct order
      pass
  elif reverse:  # Layers in reverse order
      log.warning('The atmospheric layers are in reversed order '
                  '(bottom-top).  Resorting to be from the top down.')
      if atm_in.radius is not None:
          atm_in.radius = np.flipud(atm_in.radius)
      atm_in.press = np.flipud(atm_in.press)
      atm_in.temp  = np.flipud(atm_in.temp)
      atm_in.mm    = np.flipud(atm_in.mm)
      atm_in.q     = np.flipud(atm_in.q)
      atm_in.d     = np.flipud(atm_in.d)
  else:
      log.error('The atmospheric layers are neither sorted from the '
                'bottom up, nor from the top down.')

  missing = [
      pyrat.phy.mplanet is None,
      pyrat.phy.gplanet is None,
      pyrat.phy.rplanet is None,
      atm.refpressure is None,
      ]

  if pyrat.atm.rmodelname is None and atm_in.radius is None \
     and pyrat.runmode != "opacity":
      log.error('Cannot compute radius profile.  Need to set a radius '
                'model (radmodel) or provide an input radius array in '
                'the atmospheric file.')

  if pyrat.atm.rmodelname is not None and pyrat.runmode != "opacity":
      if np.any(missing[0:3]):
          log.error('Cannot compute hydrostatic-equilibrium radius profile.  '
              'Must define at least two of mplanet, rplanet, or gplanet.')
      if missing[3]:
          log.error('Cannot compute hydrostatic-equilibrium radius profile.  '
              'Undefined reference pressure level (refpressure).')

  if pyrat.atm.rmodelname is not None and not np.any(missing):
      # Calculate hydostatic-equilibrium radius profile:
      log.msg('Reference pressure: {:.3e} {:s}.'.
              format(atm.refpressure/pt.u(atm.punits), atm.punits), indent=2)
      log.msg('Reference radius: {:8.1f} {:s}.'.
              format(pyrat.phy.rplanet/pt.u(atm.runits), atm.runits), indent=2)
      if not np.isinf(pyrat.phy.rhill):
          log.msg('Hill radius:      {:8.1f} {:s}.'.
                  format(pyrat.phy.rhill/pt.u(atm.runits), atm.runits),
                  indent=2)
      atm_in.radius = pyrat.hydro(atm_in.press, atm_in.temp, atm_in.mm,
                                  pyrat.phy.gplanet, pyrat.phy.mplanet,
                                  atm.refpressure, pyrat.phy.rplanet)

  # Check if Hydrostatic Eq. breaks down:
  ibreak = 0  # Index and flag at the same time
  if np.any(np.ediff1d(atm_in.radius) > 0.0):
      ibreak = np.where(np.ediff1d(atm_in.radius)>0)[0][0] + 1

  # Set the interpolating function (for use later):
  try:
      radinterp = sip.interp1d(atm_in.press [ibreak:],
                               atm_in.radius[ibreak:], kind='slinear')
      pressinterp = sip.interp1d(np.flipud(atm_in.radius[ibreak:]),
                             np.flipud(atm_in.press [ibreak:]), kind='slinear')
  except:
      pass

  # Set radius/pressure boundaries if exist:
  if atm.ptop is not None:
      if atm.ptop >= atm_in.press[ibreak]:
          ibreak = 0   # Turn-off break flag
  elif atm.radhigh is not None:
      if atm.radhigh <= atm_in.radius[ibreak]:
          ibreak = 0   # Turn-off break flag
          atm.ptop = pressinterp(atm.radhigh)[0]
      #else:
      #  out-of-bounds error
  else:
      atm.ptop = np.amin(atm_in.press)

  if atm.radlow is not None:
      atm.pbottom = pressinterp(atm.radlow)[0]
  elif atm.pbottom is None:
      atm.pbottom = np.amax(atm_in.press)

  if ibreak != 0 and np.isinf(pyrat.phy.rhill):
      log.error('Unbounded atmosphere.  Hydrostatic-equilibrium radius '
          'solution diverges at pressure {:.3e} bar.  Set mstar and smaxis '
          'to define a Hill radius (top boundary) and avoid error.'.
          format(atm_in.press[ibreak]/pc.bar))

  # Out of bounds errors:
  if atm.ptop < np.amin(atm_in.press) or atm.ptop > np.amax(atm_in.press):
      log.error('Top-pressure boundary (ptop={:.2e} {:s}) lies outside '
          'of the atmospheric-file range {:.2e}--{:.2e} {:s}.'.
          format(atm.ptop/pt.u(atm.punits), atm.punits,
                 np.amin(atm_in.press)/pt.u(atm.punits),
                 np.amax(atm_in.press)/pt.u(atm.punits), atm.punits))
  if atm.pbottom < np.amin(atm_in.press) or atm.pbottom > np.amax(atm_in.press):
      log.error('Bottom-pressure boundary (pbottom={:.2e} {:s}) lies outside '
          'of the atmospheric-file range {:.2e}--{:.2e} {:s}.'.
          format(atm.pbottom/pt.u(atm.punits), atm.punits,
                 np.amin(atm_in.press)/pt.u(atm.punits),
                 np.amax(atm_in.press)/pt.u(atm.punits), atm.punits))
  if atm.pbottom <= atm.ptop:
      log.error('Bottom-layer pressure ({:.2e} {:s}) must be higher than the '
                'top-layer pressure ({:.2e} {:s}).'.
                format(atm.pbottom/pt.u(atm.punits), atm.punits,
                       atm.ptop/pt.u(atm.punits), atm.punits))

  log.msg('User pressure boundaries: {:.2e}--{:.2e} bar.'.
          format(atm.ptop/pc.bar, atm.pbottom/pc.bar), indent=2)

  # Resample to equispaced log-pressure array if requested:
  if pyrat.inputs.nlayers is not None:
      atm.press = np.logspace(np.log10(atm.ptop), np.log10(atm.pbottom),
                              atm.nlayers)
      atm.radius = radinterp(atm.press)
      resample = True

  # Resample to equispaced radius array if requested:
  elif atm.radstep is not None:
      # Make equispaced radius array:
      if atm.radlow is None:
          atm.radlow  = radinterp(atm.pbottom)[0]
      if atm.radhigh is None:
          atm.radhigh = radinterp(atm.ptop)[0]

      atm.radius = np.arange(atm.radlow, atm.radhigh, atm.radstep)
      atm.nlayers = len(atm.radius)
      # Interpolate to pressure array:
      atm.press = pressinterp(atm.radius)
      resample = True

  else:  # Take the atmospheric-file sampling:
      # Get top-bottom indices:
      ilow  = np.where(atm_in.press >= atm.ptop)   [0][ 0]
      ihigh = np.where(atm_in.press <= atm.pbottom)[0][-1]
      # Take values within the boundaries:
      atm.press = atm_in.press[ilow:ihigh+1]
      atm.temp  = atm_in.temp [ilow:ihigh+1]
      atm.mm    = atm_in.mm   [ilow:ihigh+1]
      atm.q     = atm_in.q    [ilow:ihigh+1]
      atm.d     = atm_in.d    [ilow:ihigh+1]
      if atm_in.radius is not None:
          atm.radius = atm_in.radius[ilow:ihigh+1]
      atm.nlayers = len(atm.press)
      resample = False

  # Check the radii that lie within Hill radius:
  if atm.radius is not None:
      atm.rtop = pt.ifirst(atm.radius < pyrat.phy.rhill, default_ret=0)
  if atm.rtop > 0:
      log.warning('The atmospheric pressure array extends beyond '
          'the Hill radius ({:.1f} km) at pressure {:.2e} bar (layer {:d}).  '
          'Extinction beyond this layer will be neglected.'.format(
          pyrat.phy.rhill/pc.km, atm_in.press[atm.rtop]/pc.bar, atm.rtop))

  # Print radius array:
  if atm.radius is not None:
      radstr = '['+', '.join('{:9.2f}'.format(k) for k in atm.radius/pc.km)+']'
      log.msg('Radius array (km) =   {:s}'.format(radstr), indent=2, si=4)
      log.msg('Valid upper/lower radius boundaries:    {:8.1f} - {:8.1f} '
              '{:s}.'.format(atm.radius[atm.rtop]/pt.u(atm.runits),
                             atm.radius[-1]/pt.u(atm.runits),
                             atm.runits), indent=2)

  log.msg('Valid lower/higher pressure boundaries: {:.2e} - {:.2e} {:s}.'.
          format(atm.press[atm.rtop]/pt.u(atm.punits),
          atm.pbottom/pt.u(atm.punits), atm.punits), indent=2)
  log.msg('Number of valid model layers: {:d}.'.
          format(atm.nlayers-atm.rtop), indent=2)

  # Interpolate to new atm-layer sampling if necessary:
  if resample:
      tempinterp = sip.interp1d(atm_in.press, atm_in.temp, kind='slinear')
      mminterp   = sip.interp1d(atm_in.press, atm_in.mm,   kind='slinear')
      atm.temp = tempinterp(atm.press)
      atm.m    =   mminterp(atm.press)
      # Interpolate abundance profiles:
      atm.q = np.zeros((atm.nlayers, pyrat.mol.nmol))
      atm.d = np.zeros((atm.nlayers, pyrat.mol.nmol))
      for i in np.arange(pyrat.mol.nmol):
          qinterp = sip.interp1d(atm_in.press, atm_in.q[:, i], kind='slinear')
          dinterp = sip.interp1d(atm_in.press, atm_in.d[:, i], kind='slinear')
          atm.q[:,i] = qinterp(atm.press)
          atm.d[:,i] = dinterp(atm.press)

  # Interpolate isotopes partition function:
  log.msg('Number of isotopes: {:d}'.format(pyrat.iso.niso), indent=2)
  # Initialize the partition-function array for pyrat.iso:
  pyrat.iso.z = np.zeros((pyrat.iso.niso, atm.nlayers))
  for db in pyrat.lt.db:            # For each Database
      for j in np.arange(db.niso):  # For each isotope in DB
          log.debug('Interpolating (isotope ID {:2d}) partition function.'.
                  format(db.iiso+j), indent=4)
          zinterp = sip.interp1d(db.temp, db.z[j], kind='slinear')
          pyrat.iso.z[db.iiso+j] = zinterp(atm.temp)

  # Base abundance profiles:
  atm.qbase = np.copy(atm.q)

  log.head('Make atmosphere done.')
