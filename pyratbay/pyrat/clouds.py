# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np


def absorption(pyrat):
  """
  Evaluate the total cloud absorption in the atmosphere.
  """
  pyrat.cloud.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for model in pyrat.cloud.models:
      if model.name == 'deck':
          model.extinction(pyrat.spec.wn, pyrat.atm.press, pyrat.atm.radius)
          pyrat.cloud.ec += model.ec
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
          model.extinction(pyrat.spec.wn, pyrat.atm.press, pyrat.atm.radius)
          ec.append(model.ec[layer])
      else:
          imol = np.where(pyrat.mol.name == model.mol)[0][0]
          ec.append(model.ec[layer] * pyrat.atm.d[layer,imol])
      label.append(model.name)
  return ec, label

