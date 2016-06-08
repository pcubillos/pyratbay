import os, re, sys
from numpy import get_include
from setuptools import setup, Extension

topdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(topdir + "/pyratbay")
import VERSION as ver

srcdir = topdir + '/src_c/'          # C-code source folder
incdir = topdir + '/src_c/include/'  # Include filder with header files

# Get all file from source dir:
files = os.listdir(srcdir)

# Filter the results for just the c files:
files = list(filter(lambda x:     re.search('.+[.]c$',     x), files))
files = list(filter(lambda x: not re.search('[.#].+[.]c$', x), files))

inc = [get_include(), incdir]
eca = ['-ffast-math']
ela = []

extensions = []
for i in range(len(files)):
  e = Extension(os.path.splitext(files[i])[0],
                sources=['{:s}{:s}'.format(srcdir, files[i])],
                include_dirs=inc,
                extra_compile_args=eca,
                extra_link_args=ela)
  extensions.append(e)


setup(name         = "pyratbay",
      version      = "{:d}.{:d}.{:d}".format(ver.PBAY_VER, ver.PBAY_MIN,
                                             ver.PBAY_REV),
      author       = "Patricio Cubillos",
      author_email = "patricio.cubillos@oeaw.ac.at",
      url          = "https://github.com/pcubillos/pyratbay",
      packages     = ["pyratbay"],
      license      = ["FINDME"],
      description  = "Python Radiative Transfer in a Bayesian Framework.",
      include_dirs = inc,
      ext_modules  = extensions)
