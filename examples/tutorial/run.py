# Preamble
# --------
# To correctly execute this script, set the correct path to the source
# code.   The path is given as if the Python session runs from a
# 'run_tutorial/' folder at the same level than the repository, i.e.:
#    parentdir/
#    |-- pyratbay/
#    `-- run_tutorial/
#  Alternatively, set the appropriate path in sys.path.append().

"""
Once you are in the run_tutorial/ folder, download the following Kurucz
stellar model file with this shell command for example (or you can
also use the browser):

curl http://kurucz.harvard.edu/grids/gridp00odfnew/fp00k2odfnew.pck -o kurucz_fp00k2odfnew.pck
"""

import sys, os
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# Edit the path to the Pyrat-Bay package if necessary:
sys.path.append("../pyratbay")
import pyratbay as pb


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Make a TLI file containing the opacity line-transition info:
"""
# NOTE: Before running the Python script, you need to download and
# unzip this HITRAN file into the working directory.
# Download from the shell prompt using wget:
wget --user=HITRAN --password=getdata -N https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip

# Alternatively, download from the shell prompt using curl:
curl -u HITRAN:getdata https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip -o 01_hit12.zip

# Unzip the file:
unzip 01_hit12.zip
"""

pb.pbay.run("tutorial_tli.cfg")


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Make temperature profiles:

# Generate an isothermal PT profile (output values in CGS units):
pressure, T_isothermal = pb.pbay.run("tutorial_pt-isothermal.cfg")
# Generate a TCEA PT profile:
pressure, T_tcea       = pb.pbay.run("tutorial_pt-tcea.cfg")

# Plot the PT profiles:
plt.figure(-1)
plt.clf()
plt.semilogy(T_isothermal, pressure/pb.constants.bar, color="b",
             lw=2, label='Isothermal')
plt.semilogy(T_tcea, pressure/pb.constants.bar, color="r",
             lw=2, label='TCEA')
plt.ylim(100, 1e-5)
plt.legend(loc="best")
plt.xlabel("Temperature  (K)")
plt.ylabel("Pressure  (bar)")
plt.savefig("pyrat_PT_tutorial.pdf")


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Make Atmospheric file:
# (Might need to install sympy (only the first time!)
#  pip2 install sympy==0.7.1 )

# Generate a TEA atmospheric model:
pressure, temperature, abundances = pb.pbay.run("tutorial_atmosphere-tea.cfg")
# Generate a uniform-abundance atmospheric model:
pressure, temperature, abundances = pb.pbay.run("tutorial_atmosphere-uniform.cfg")

# Read the atmospheric files:
spec, press, temp, q_tea     = pb.atmosphere.readatm("WASP-00b.atm")
spec, press, temp, q_uniform = pb.atmosphere.readatm("WASP-00c.atm")

# Plot the results:
plt.figure(-2)
plt.clf()
ax = plt.subplot(211)
for i in np.arange(len(spec)):
  plt.loglog(q_tea[:,i], press, label=spec[i], lw=2)

plt.ylim(np.amax(press), np.amin(press))
plt.xlim(1e-10, 1.0)
plt.legend(loc='best', fontsize=11)
plt.ylabel("Pressure  (bar)")
ax = plt.subplot(212)
for i in np.arange(len(spec)):
  plt.loglog(q_uniform[:,i], press, label=spec[i], lw=2)

plt.ylim(np.amax(press), np.amin(press))
plt.xlim(1e-10, 1.0)
plt.xlabel("Mole mixing fraction")
plt.ylabel("Pressure  (bar)")
plt.savefig("pyrat_atmosphere_tutorial.pdf")


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Compute spectrum:

pyrat = pb.pbay.run("tutorial_spectrum.cfg")

# Plot:
plt.figure(-3)
plt.clf()
ax = plt.subplot(111)
plt.semilogx(1e4/pyrat.spec.wn, pyrat.spec.spectrum, "b-")
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0])
plt.xlim(0.3, 5.0)
plt.ylabel("Modulation  (Rp/Rs)^2")
plt.xlabel("Wavelength  (um)")

# Alternative plotting:
ax = pb.plots.spectrum(pyrat=pyrat)
ax.set_xscale('log')
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0])
plt.show()
plt.savefig("pyrat_transmission-spectrum_tutorial.pdf")

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generate opacity grid:

pyrat = pb.pbay.run("tutorial_opacity.cfg")


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# A simple MCMC:

pyrat = pb.pbay.run("tutorial_mcmc.cfg")
