from numpy import get_include
import os, re, sys
from distutils.core import setup, Extension

srcdir = './'    # C-code source folder
libdir = 'lib/'  # Where the shared objects are put

# Get all file from source dir:
files = os.listdir(srcdir)
# Filter the results for just the c files:
files = filter(lambda x:     re.search('.+[.]c$',     x), files)
files = filter(lambda x: not re.search('[.#].+[.]c$', x), files)

ext_mod = []
inc = [get_include()]
ext_comp_args = ", ".join([#"'-fopenmp'",
                           "'-ffast-math'"
                                          ])
ext_link_args = ", ".join([#"'-lgomp'"
                                     ])

for i in range(len(files)):
  exec("mod{:d} = Extension('{:s}', "
                           "sources=['{:s}{:s}'], "
                           "include_dirs=inc, "
                           "extra_compile_args=[{:s}], "
                           "extra_link_args=[{:s}])".
        format(i, files[i].rstrip('.c'), srcdir, files[i],
               ext_comp_args, ext_link_args))

  exec('ext_mod.append(mod{:d})'.format(i))

setup(name=libdir, version='1.0', description='c extension functions',
      ext_modules = ext_mod)
