# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import re
import sys
import setuptools
from setuptools import setup, Extension

from numpy import get_include

sys.path.append(os.path.join(os.path.dirname(__file__), 'pyratbay'))
from VERSION import __version__


srcdir = 'src_c/'          # C-code source folder
incdir = 'src_c/include/'  # Include filder with header files

cfiles = os.listdir(srcdir)
cfiles = list(filter(lambda x: re.search('.+[.]c$', x), cfiles))
cfiles = list(filter(lambda x: not re.search('[.#].+[.]c$', x), cfiles))

inc = [get_include(), incdir]
eca = ['-ffast-math']
ela = []

extensions = [
    Extension(
        'pyratbay.lib.' + cfile.rstrip('.c'),
        sources=[f'{srcdir}{cfile}'],
        include_dirs=inc,
        extra_compile_args=eca,
        extra_link_args=ela)
    for cfile in cfiles
    ]

long_description = """
.. image:: https://raw.githubusercontent.com/pcubillos/pyratbay/master/docs/figures/pyrat_logo.png
   :width: 50%

|Build Status|  |docs|  |PyPI|  |conda|  |License|

``Pyrat Bay``: Python Radiative Transfer in a Bayesian framework


Install as:

.. code-block:: shell

    pip install pyratbay

Or:

.. code-block:: shell

    conda install -c conda-forge pyratbay

Docs at: https://github.com/pcubillos/pyratbay

Cite as:

.. code-block:: bibtex

  @ARTICLE{CubillosBlecic2021mnrasPyratBay,
         author = {{Cubillos}, Patricio E. and {Blecic}, Jasmina},
          title = "The {Pyrat Bay} Framework for Exoplanet Atmospheric Modeling: A Population Study of Hubble/WFC3 Transmission Spectra",
           year = 2021,
        journal = {\mnras},
            doi = {10.1093/mnras/stx0000},
         adsurl = {https://ui.adsabs.harvard.edu/abs/2021MNRAS.000.0000C},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System},
  }
.. |Build Status| image:: https://travis-ci.com/pcubillos/pyratbay.svg?branch=master
   :target: https://travis-ci.com/pcubillos/pyratbay

.. |docs| image:: https://readthedocs.org/projects/pyratbay/badge/?version=latest
    :target: https://pyratbay.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |PyPI| image:: https://img.shields.io/pypi/v/pyratbay.svg
    :target:      https://pypi.org/project/pyratbay/
    :alt: Latest Version

.. |conda| image:: https://img.shields.io/conda/vn/conda-forge/pyratbay.svg
    :target: https://anaconda.org/conda-forge/pyratbay

.. |License| image:: https://img.shields.io/github/license/pcubillos/pyratbay.svg?color=blue
    :target: https://pyratbay.readthedocs.io/en/latest/license.html
"""


package_data = {'pyratbay': [
    'data/AsplundEtal2009.txt',
    'data/atoms.dat',
    'data/molecules.dat',
    'data/isotopes.dat',
    'data/TEA_gdata_defaults.txt',
    'data/tips_2017.pkl',
    'data/filters/*.dat',
    'data/CIA/*.dat',
    'TEA/lib/abundances.txt',
    'TEA/lib/stoich.txt',
    'TEA/lib/TEA.cfg',
    'TEA/lib/gdata/*.txt',
    'TEA/*/*.py',
    ]}


setup(
    name = 'pyratbay',
    version = __version__,
    author = 'Patricio Cubillos',
    author_email = 'patricio.cubillos@oeaw.ac.at',
    url = 'https://github.com/pcubillos/pyratbay',
    packages = setuptools.find_packages(),
    package_data = package_data,
    install_requires = [
        'numpy>=1.8.1',
        'scipy>=0.13.3',
        'matplotlib>=1.3.1',
        'sympy>=0.7.6',
        'mc3>=3.0.7',
        ],
    tests_require = [
        'pytest>=3.9',
        'scipy>=1.4.1',
        ],
    license = 'GPLv2',
    description = 'Python Radiative Transfer in a Bayesian Framework.',
    long_description = long_description,
    long_description_content_type = 'text/x-rst',
    include_dirs = inc,
    entry_points = {'console_scripts': ['pbay = pyratbay.__main__:main']},
    ext_modules = extensions,
    )
