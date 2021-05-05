# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import numpy as np


def absorption(pyrat):
  """
  Evaluate the total cloud absorption in the atmosphere.
  """
  pyrat.cloud.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for model in pyrat.cloud.models:
      if model.name == 'deck':
          model.extinction(pyrat.atm.press, pyrat.atm.radius, pyrat.atm.temp)
          continue

      # Calculate the extinction coefficient (in cm2 molecule-1):
      model.extinction(pyrat.spec.wn, pyrat.atm.press)
      imol = np.where(pyrat.mol.name == model.mol)[0][0]
      # Densities in molecules cm-3:
      dens = pyrat.atm.d[:,imol]
      # Cloud absorption (cm-1):
      pyrat.cloud.ec += model.ec * np.expand_dims(dens, axis=1)


def get_ec(pyrat, layer):
  """
  Extract per-model extinction coefficient at requested layer.
  """
  ec, label = [], []
  for model in pyrat.cloud.models:
      if model.name == 'deck':
          model.extinction(pyrat.atm.press, pyrat.atm.radius, pyrat.atm.temp)
          # Note this is just a filler, the actual code does not use this:
          e = np.zeros(pyrat.spec.nwave) + int(layer > model.itop)
          ec.append(e)
      else:
          imol = np.where(pyrat.mol.name == model.mol)[0][0]
          ec.append(model.ec[layer] * pyrat.atm.d[layer,imol])
      label.append(model.name)
  return ec, label

