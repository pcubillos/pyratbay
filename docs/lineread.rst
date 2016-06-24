.. _lineread:

Lineread Module
===============

The ``lineread`` package formats online-available line-by-line opacity
data into a line-transition information (TLI) file, which is used as
input for ``pyrat``.


Quick Example
-------------

Pyrat Bay (version 1.X) is known to work (at least) on Unix/Linux (Ubuntu)
and OSX (10.9+) machines, with the following software:

* Python (version 2.7)
* Numpy (version 1.8.2+)
* Scipy (version 0.13.3+)
* Matplotlib (version 1.3.1+)

MC3 may work with previous versions of these software.
However we do not guarantee nor provide support for that.

Walkthrough
-----------

To obtain the latest stable Pyrat-Bay code, clone the repository to
your local machine with the following terminal commands.
First, from any given directory, clone the Pyrat-Bay repository:

.. code-block:: shell

  topdir=`pwd`
  git clone https://github.com/pcubillos/pyratbay

Quick Example: ``pyrat`` forward modeling
-----------------------------------------

The following script quickly lets you calculate a water transmission
spectrum between 1 and 5.5 um.  These instructions are meant to be
executed from a Shell terminal.  To begin, follow the instructions
in the previous Section to install and compile the code.
Now, create a working directory to place the files and execute the programs:

.. code-block:: shell

   cd $topdir
   mkdir run01
   cd run01

Download the water line-transition database from the HITRAN server with wget or curl, and unpack the zip file:

.. code-block:: shell

   # Using wget:
   wget --user=HITRAN --password=getdata -N https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip
   # Or alternatively, curl:
   curl -u HITRAN:getdata https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip -o 01_hit12.zip
   unzip 01_hit12.zip

Copy the Lineread configuration file and run Lineread to make the Transition-Line-Information (TLI) file:

.. code-block:: shell

   cp $topdir/pyratbay/examples/lineread_demo/demo_hitran.cfg .
   $topdir/pyratbay/lineread.py -c demo_hitran.cfg

Copy the Pyrat configuration file and run it to compute the spectrum:

.. code-block:: shell

   cp $topdir/Pyrat-Bay/examples/pyrat_demo/* .
   $topdir/Pyrat-Bay/pyrat.py -c demo_transit.cfg

Outputs
^^^^^^^

That's it, now let's see the results.  Run this Python script:

.. code-block:: python

  import matplotlib.pyplot as plt
  import sys
  sys.path.append("../Pyrat-Bay/scripts/")
  import readpyrat as rp
  wlength, flux = rp.readspectrum("transit_spectrum.dat", 0)
  
  plt.figure(0, (8,5))
  plt.clf()
  plt.title("Water Modulation Spectrum")
  plt.plot(wlength, flux, "b", label="H2O/H2/He Model")
  plt.xlabel("Wavelength  (um)")
  plt.ylabel("Modulation")
  plt.show()


Pyrat-Bay will print out to screen some stuff:

.. code-block:: none

   Start MCMC chains  (Tue Jan  5 13:11:22 2016)
   
   ...
  
   [::        ]  20.0% completed  (Tue Jan  5 13:11:22 2016)
   Out-of-bound Trials:
    [0 0 0]
   Best Parameters:   (chisq=87.5664)
   [ 2.81119952 -2.33026943  0.48622898]

   ...
