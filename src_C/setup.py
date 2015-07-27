from numpy import get_include
import os, re, sys
from distutils.core import setup, Extension

srcdir = './'           # C-code source folder
incdir = '../include/'  # Include filder with header files

# Get all file from source dir:
files = os.listdir(srcdir)

# Filter the results for just the c files:
files = list(filter(lambda x:     re.search('.+[.]c$',     x), files))
files = list(filter(lambda x: not re.search('[.#].+[.]c$', x), files))

inc = [get_include(), incdir]
eca = ['-ffast-math']  # '-fopenmp'
ela = []               # '-lgomp'

extensions = []
for i in range(len(files)):
  e = Extension(files[i].rstrip('.c'),
                sources=['{:s}{:s}'.format(srcdir, files[i])],
                include_dirs=inc,
                extra_compile_args=eca,
                extra_link_args=ela)
  extensions.append(e)


setup(name         = "Pyrat-Bay C",
      version      = '1.0',
      author       = "Patricio Cubillos",
      author_email = "pcubillos@fulbrightmail.org",
      url          = "https://github.com/pcubillos/Pyrat-Bay",
      description  = "Pyrat-Bay C extension function",
      ext_modules  = extensions)
