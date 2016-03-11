import sys, os, shutil
import subprocess
import scipy.constants as sc

pbpath = "../Pyrat-Bay"

sys.path.append(pbpath)
import pyratbay.pyratbay as pb
import pyratbay.pyrat    as p

sys.path.append(pbpath+"/lib")
import pt as pt

sys.path.append(pbpath+"/scripts")
import readatm as ra

# Directory of BART.py file:
#PBdir  = os.path.dirname(os.path.realpath(__file__))
TEAdir = pbpath + "/modules/TEA/"

# Pressure array (bars):
pbottom = 100.0
ptop    = 1e-8
nlayers = 150
pressure = np.logspace(np.log10(ptop), np.log10(pbottom), nlayers) 

# WASP-126 parameters:
Rstar   = 1.27*6.955e8 # m
Tstar   = 5800.0       # K
Tint    =  100.0       # K
gplanet =  676.1       # cm s-2
smaxis  = 0.045*sc.au  # m

# Fitting parameters: [log10(kappa), log10(g1), log10(g2), alpha, beta]
params = np.asarray(  [-1.0,         -0.25,     -0.8,      0.0,   1.0])

# Calculate the temperature profile:
temp = pt.TCEA(params, pressure, Rstar, Tstar, Tint, smaxis, gplanet)

# Generate elemental-abundances file:
solar  = pbpath + "/inputs/AsplundEtal2009.txt"
afile  = "atomic.atm"
xsolar = 1.0
swap   = None
pb.ma.makeatomic(solar, afile, xsolar, swap)


# Calculate the temperature profile:
patm = "preatm.tea"
# Input elemental (atomic) composition:
elements = ["H", "He", "C", "N", "O", "Na", "K"]
# Output species composition:
species = ["H", "He", "C", "N", "O", "Na", "K", "H2",
           "CO", "CO2", "CH4", "H2O", "NH3", "C2H2", "C2H4", "HCN"]
pb.ma.makepreatm(pressure, temp, afile, elements, species, patm)


# Make a TEA config file:
pb.mc.makeTEA(abun_file=afile)
# Execute TEA:
proc = subprocess.Popen([TEAdir + "tea/runatm.py", patm, 'TEA'])
# mv TEA.atm file
atmfile = "WASP-126b.atm"
pb.ma.TEA2pyrat("TEA/TEA/results/TEA.tea", atmfile)
shutil.rmtree("TEA")

# See the results:
spec, press, temp, q = ra.readatm(atmfile)
fmt = ["b", "c", "y", "m", "k", "r--", "g--", "r", "k:", "c:", "m:", "b:", "y:", "r:", "g:", "b+"]

plt.figure(0)
plt.clf()
for i in np.arange(len(spec)):
  plt.loglog(q[:,i], pressure, fmt[i], label=spec[i])
  plt.ylim(press[-1], press[0])

plt.legend(loc='best')

