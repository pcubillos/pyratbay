Pyrat-Bay
=========

> Python Radiative Transfer in a Bayesian Framework.

### Table of Contents:
* [Team Members](#team-members)
* [Getting Started](#getting-started)
* [Install and Compile](#install-and-compile)
* [Quick Example](#quick-example)
* [Be Kind](#be-kind)
* [License](#license)

### Team Members:
* [Patricio Cubillos](https://github.com/pcubillos/) (Space Research Institute, Graz) <pcubillos@fulbrightmail.org>

### Getting Started:
Take a look at the Pyrat-Bay [Wiki](https://github.com/pcubillos/Pyrat-Bay/wiki).

### Install and Compile:
To obtain the Pyrat-Bay code download the latest stable version from the releases page (TBD). Alternatively, clone the repository to your local machine with the following terminal commands.  First, create a top-level directory to place the code:  
```shell
mkdir pyrat_demo/  
cd pyrat_demo/  
topdir=`pwd`
```

Clone the repository to your working directory:  
```shell
git clone --recursive  https://github.com/pcubillos/Pyrat-Bay
```

Compile the C programs:  
```shell
cd $topdir/Pyrat-Bay/src_C/
make
cd $topdir/Pyrat-Bay/ctips/
make
```

To remove the program binaries, execute (from the respective directories):  
```shell
make clean
```

### Quick Example:
The following script quickly lets you calculate a water transmission spectrum between 1 and 5.5 um.  These instructions are meant to be executed from a Shell terminal.  To begin, follow the instructions in the previous Section to install and compile the code.  Now, create a working directory to place the files and execute the programs:
```shell
cd $topdir
mkdir run/  
cd run/  
```

Download the water line-transition database from the HITRAN server with wget:
```shell
wget --user=HITRAN --password=getdata -N https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip
unzip 01_hit12.zip
```
Or alternatively, use curl:
```shell
curl -u HITRAN:getdata https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip -o 01_hit12.zip
unzip 01_hit12.zip
```

Copy the Lineread configuration file and run Lineread to make the Transition-Line-Information (TLI) file:
```shell
cp $topdir/Pyrat-Bay/examples/lineread_demo/demo_hitran.cfg .
$topdir/Pyrat-Bay/lineread.py -c demo_hitran.cfg
```

Copy the Pyrat configuration file and run it to compute the spectrum:
```shell
cp $topdir/Pyrat-Bay/examples/pyrat_demo/* .
$topdir/Pyrat-Bay/pyrat.py -c demo_transit.cfg
```

To check out the results, run this Python script:
```python
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
```

### Be Kind:
Please reference these papers if you found this module useful for your research:  
  [Cubillos et al. 2016: Yet Another Open-source Radiative-Transifer Code for Exoplanet Modeling](), in preparation.   
Thanks!


### License:

There should be a license soon, right?
