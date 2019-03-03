import os
import re
import sys
import setuptools
from setuptools import setup, Extension

from numpy import get_include

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
for efile in files:
    #e = Extension(os.path.splitext(efile)[0],
    e = Extension('pyratbay.lib.'+efile.rstrip('.c'),
                  sources=['{:s}{:s}'.format(srcdir, efile)],
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
      packages     = setuptools.find_packages(),
      install_requires = ['numpy>=1.8.1',
                          'scipy>=0.13.3',
                          'matplotlib>=1.3.1',
                          'sympy>=0.7.6'],
      license      = "TBD",
      description  = "Python Radiative Transfer in a Bayesian Framework.",
      include_dirs = inc,
      ext_modules  = extensions)
