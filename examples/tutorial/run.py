# Preamble
# --------
# To correctly execute this script, set the correct path to the Pyrat Bay
# source code.  The path is given as if the Python session runs from a
# 'run/' folder at the same level than the repository folder, as in:
#    rootdir/
#    |-- pyratbay/
#    `-- run/
#  Alternatively, set pbpath to the appropriate path.

import sys, os
import numpy as np
import matplotlib.pyplot as plt
# This line allows you to keep a plotting window open:
plt.ion()

# Set the path to the pyratbay dir:
pbpath = "../pyratbay"

# Import the Pyrat Bay package:
sys.path.append(pbpath)
import pyratbay as pb

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Make a TLI file with opacity line-transition info:
pb.pbay.run("tutorial_lineread.cfg")

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Make Temperature profile: 

# Generate PT profile (output values in CGS units):
pressure, temperature = pb.pbay.run("tutorial_pt.cfg")

# Plot the resulting profile:
plt.figure(0)
plt.clf()
plt.semilogy(temperature, pressure/pb.constants.bar, color="b", lw=2)
plt.ylim(100, 1e-5)
plt.xlabel("Temperature  (K)")
plt.ylabel("Pressure  (bar)")
plt.show()

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Make Atmospheric file:
# (Might need to install sympy (only the first time!)
#  pip2 install sympy==0.7.1 )

# Generate a TEA atmospheric model:
pressure, temperature, abundances = pb.run("tutorial_atmosphere.cfg")

# Read the atmospheric file:
spec, press, temp, q = pb.atmosphere.readatm("WASP-00b.atm")

# Plot the results:
plt.figure(1)
plt.clf()
for i in np.arange(len(spec)):
  plt.loglog(q[:,i], press, label=spec[i], lw=2)

plt.ylim(np.amax(press), np.amin(press))
plt.xlim(1e-10, 1.0)
plt.legend(loc='best')
plt.xlabel("Mole mixing fraction")
plt.ylabel("Pressure  (bar)")

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Compute spectrum:

pyrat = pb.pbay.run("tutorial_spectrum.cfg")

plt.figure(2)
plt.clf()
ax = plt.subplot(111)
plt.semilogx(1e4/pyrat.spec.wn, pyrat.spec.spectrum, "b-")
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0])
plt.xlim(0.3, 5.0)
plt.ylabel("Modulation spectrum  (Rp/Rs)^2")
plt.xlabel("Wavelength  (um)")
plt.show()

ax = pb.plots.spectrum(pyrat=pyrat)
ax.set_xscale('log')
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0])
plt.show()



