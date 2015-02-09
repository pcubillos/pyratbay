# ****************************** START LICENSE ******************************
# ******************************* END LICENSE ******************************

import sys, os
import numpy as np
import subprocess as sp

import ptools     as pt
import pconstants as pc
from driver import dbdriver

# Directory of db_hitran.py:
DBHdir = os.path.dirname(os.path.realpath(__file__))

class hitran(dbdriver):
  def __init__(self, dbfile, pffile):
    """
    Initialize the basic database info.

    Parameters:
    -----------
    dbfile: String
       File with the Database info as given from HITRAN.
    pffile: String
       File with the partition function.

    Modification History:
    ---------------------
    2014-08-01  patricio  Added documentation.
    """
    super(hitran, self).__init__(dbfile, pffile)

    self.recsize   = 162 # Record length
    self.recwnpos  =   3 # Wavenumber     position in record
    self.recisopos =   2 # Isotope        position in record
    self.reclinpos =  15 # Line intensity position in record
    self.recApos   =  25 # Einstein coef  position in record
    self.recairpos =  35 # Air broadening position in record
    self.recelpos  =  45 # Low Energy     position in record
    self.recg2pos  = 155 # Low stat weight position in record
    self.recmollen =   2 # Molecule   record length
    self.recwnlen  =  12 # Wavenumber record length
    self.reclinend =  25 # Line intensity end position
    self.recelend  =  55 # Low Energy     end position
    self.T0       = 296.0 # K

    self.molID = self.getMolec()
    # Get info from HITRAN configuration file:
    self.molecule, self.isotopes, self.mass, self.isoratio, self.gi = \
                                                    self.getHITinfo(self.molID)
    # Database name: 
    self.name = "HITRAN " + self.molecule


  def readwl(self, dbfile, irec):
    """
    Read wavelength parameter from irec record in dbfile database.

    Parameters:
    -----------
    dbfile: File object
       File where to extract the wavelength.
    irec: Integer
       Index of record.
    dbtype: String
       Database type: ['ps', 'hit']

    Returns:
    --------
    rec_wl: Float
       Wavelength value at record irec, as given in dbfile database.

    Modification History:
    ---------------------
    2014-03-05  patricio  Initial implementation, based on Madison's
                          code.          pcubillos@fulbrightmail.org
    2014-03-10  patricio  Updated dbtype to match command-line-argument
                          sytax.  Updated HITRAN data type.
    """
    # Set pointer at required wavenumber record:
    dbfile.seek(irec*self.recsize + self.recwnpos)
    # Read:
    wavenumber = dbfile.read(self.recwnlen)
    # Convert to float:
    rec_wl = float(wavenumber)

    return rec_wl


  def getpf(self, verbose):
    """                                    
    Calculate the partition function for a grid of temperatures for
    HITRAN molecules.

    Notes:
    ------
    This function is a wrapper of the Fortran TIPS routine written
    by R. R. Gamache.
    The range of temperatures is limited to: 70K -- 3000K.

    Parameters:
    -----------
    dbfile: String
       Database filename.
    pffile: String
       Partition Function filename.
    verbose: Integer
       Verbosity threshold.

    Returns:
    --------
    PF: 2D ndarray
       A 2D array (N temperatures, N isotopes + 1) with the partition
      function values for each isotope (columns) as function of temperature
       (first column).
    DBname: String
       The database name.
    isoNames: List
       List with the P&S isotope names.
    mass: List
       List with the P&S isotope masses.

    Modification History:
    ---------------------
    2012-11-22  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    2014-03-10  patricio  Adapted previous code.
    """

    # Get molecule ID as  integer:
    molID = int(self.molID)
 
    # Get Number of isotopes:
    Niso = len(self.isotopes)

    # Array of temperatures:
    Temp = np.arange(70.0, 3000.1, 10)
    # TIPS range of temperature is: [70K, 3000K]
    Ntemp = len(Temp)

    # Output array for table of Temperature and PF values: 
    PF = np.zeros((Niso, Ntemp), np.double)

    for i in np.arange(Niso):
      # Molecule ID, isotope ID, and temperature array as strings:
      smol  = str(molID) + "\n"
      siso  = str(self.isotopes[i]) + "\n"
      stemp = np.asarray(Temp, '|S6')

      # String input for TIPS:
      input = smol + siso + ("\n" + smol + siso).join(stemp) + "\n99\n"

      # Call fortran routine to calculate the partition function:
      p = sp.Popen([self.pffile], stdout=sp.PIPE, stdin=sp.PIPE,
                   stderr=sp.STDOUT)
      #print("Size: %d"%np.size(PF[i]))
      #tmp =  p.communicate(input=input)[0].strip().split("\n")
      #print(tmp)
      PF[i] = p.communicate(input=input)[0].strip().split("\n") 
      PF[i] /= self.gi[i]
    return Temp, PF


  def getMolec(self):
    """
    Get the HITRAN molecule index

    Parameters:
    -----------
    dbfile: String
       The data filename

    Modification History:
    ---------------------
    2014-03-21  patricio  Initial implementation.
    """
    # Open file and read first two characters:
    data = open(self.dbfile, "r")
    molID  = data.read(self.recmollen)
    data.close()
    # Set database name:
    return molID #self.molname[molID-1]


  def getHITinfo(self, molID):
    """
    Get HITRAN info from configuration file.

    Parameters:
    -----------
    molID: Integer
       Molecule ID as given in HITRAN database

    Returns:
    --------
    molname:  Molecule's name
    isotopes: Isotopes names
    mass:     Isotopes mass
    isoratio: Isotopic abundance ratio
    gi:       State-independent statistical weight

    """
    # Read HITRAN configuration file from inputs folder:
    hfile = open(DBHdir + '/../inputs/hitran.dat', 'r')
    lines = hfile.readlines()
    hfile.close()
  
    isotopes = []
    mass     = []
    isoratio = []
    gi       = []
 
    # Get values for our molecule: 
    for i in np.arange(len(lines)):
      if lines[i][0:2] == molID:
        line = lines[i].split()
        molname  = line[1]
        gi.      append(  int(line[3]))
        isotopes.append(      line[2] )
        isoratio.append(float(line[4]))
        mass.    append(float(line[5]))

    return molname, isotopes, mass, isoratio, gi    

 
  def PFzero(self):
    """
    Calculate the partition function for the isotopes at T0 = 296K.

    Notes:
    ------
    This function is a wrapper of the Fortran TIPS routine written
    by R. R. Gamache.
    The range of temperatures is limited to: 70K -- 3000K.

    Parameters:
    -----------
    pffile: String
       Partition Function filename.
    molID: Integer
       Molecule ID as given in HITRAN database.

    Returns:
    --------
    PFzero: 1D ndarray
       An array (N isotopes) with the partition function evaluated at
       T0 = 296 K for each isotope.

    Modification History:
    ---------------------
    2012-03-10  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    """
    # Get number of isotopes:
    Niso = len(self.isotopes)

    # Output array for table of Temperature and PF values:
    PFzero = np.zeros(Niso, np.double)

    # Molecule ID, isotope ID, and temperature array as strings:
    smol  = self.molID + "\n"
    siso  = np.asarray(self.isotopes, '|S5')
    stemp = "\n" + str(self.T0) + "\n"

    # String input for TIPS:
    input = smol + ("\n" + stemp + smol).join(siso) + stemp + "\n99\n"

    # Call fortran routine to calculate the partition function:
    p = sp.Popen([self.pffile], stdout=sp.PIPE, stdin=sp.PIPE, stderr=sp.STDOUT)
    PFzero[:] = p.communicate(input=input)[0].strip().split("\n")

    return PFzero / self.gi


  def dbread(self, iwn, fwn, verbose, *args):
    """
    Read a HITRAN or HITEMP database (dbfile) between wavenumbers iwn and fwn.

    Parameters:
    -----------
    dbfile: String
       A HITRAN or HITEMP database filename.
    iwl: Scalar
       Initial wavenumber limit (in cm-1).
    fwl: Scalar
       Final wavenumber limit (in cm-1).
    verbose: Integer
       Verbosity threshold.
    pffile: String
       Partition function filename.

    Returns:
    --------
    wnumber: 1D ndarray (double)
      Line-transition central wavenumber (cm-1).
    gf: 1D ndarray (double)
      gf value (unitless).
    elow: 1D ndarray (double)
      Lower-state energy (cm-1).
    isoID: 2D ndarray (integer)
      Isotope index (1, 2, 3, ...).

    Modification History:
    ---------------------
    2012-12-10  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    2014-03-10  patricio  Adapted for pylineread.
    2014-07-06  patricio  Updated to return 1D arrays.
    2015-02-01  patricio  Changed wavelengths to wavenumbers.
    """
    # Get Total number of transitions in file:
    data = open(self.dbfile, "r")
    data.seek(0, 2)
    nlines   = data.tell() / self.recsize
    # Get Molecule ID:
    molID = int(self.molID)

    # Get database limiting wavenumbers:
    lastwn  = self.readwl(data, nlines-1)
    firstwn = self.readwl(data, 0)

    # Find the record index for fwn:
    if fwn < lastwn:
      irec_init = self.binsearch(data, fwn, 0, nlines,    1)
    else:
      irec_init = nlines-1
    # Find the record index for iwn:
    if iwn > firstwn:
      irec_fin  = self.binsearch(data, iwn, 0, irec_init, 0)
    else:
      irec_fin  = 0

    # Allocate arrays for values to extract:
    wnumber = np.zeros(irec_init-irec_fin+1, np.double)
    gf      = np.zeros(irec_init-irec_fin+1, np.double)
    elow    = np.zeros(irec_init-irec_fin+1, np.double)
    isoID   = np.zeros(irec_init-irec_fin+1,       int)
    # Einstein A coefficient:
    A21     = np.zeros(irec_init-irec_fin+1, np.double)
    # Lower statistical weight:
    g2      = np.zeros(irec_init-irec_fin+1, np.double)
    # Line intensity, used to calculate gf:
    S0      = np.zeros(irec_init-irec_fin+1, np.double)
    # Wavelength:
    wlength = np.zeros(irec_init-irec_fin+1, np.double)

    # Get the partition function at T0:
    PFzero = self.PFzero()

    i = 0  # Stored record index
    chk = 0
    interval = float((irec_fin - irec_init)/20)  # Check-point interval
    while irec_init - i >= irec_fin:
      data.seek( (irec_init-i) * self.recsize )
      # Read in wavenumber
      line = data.read(self.recsize)
      # Get the isotope index:
      isoID[i]   = float(line[self.recisopos:self.recwnpos ])
      # Get the wavenumber:
      wnumber[i] = float(line[self.recwnpos: self.reclinpos])
      # Get the line intensity:
      S0[i]      = float(line[self.reclinpos:self.reclinend])
      # Get Elow:
      elow[i]    = float(line[self.recelpos: self.recelend ])
      A21[i]     = float(line[self.recApos:  self.recairpos])
      g2[i]      = float(line[self.recg2pos: self.recsize  ])
      # Print a checkpoint statement every 1/20th interval
      if verbose > 1:
        pos = float(data.tell()/self.recsize)
        if (pos/interval)%1 == 0.0:
          pt.msg(verbose-1, "checkpoint {:d}/20...".format(chk), 2)
          chk += 1
          pt.msg(verbose-3, "Wavenumber: {:.2f}, S0: {:.3e}, Elow: {:.3e}".
                            format(wnumber[i], S0[i], elow[i]), 2)
          pt.msg(verbose-3, "Wavelength: {:.3f}, IsoID: {:d}, Elow: {:.3f}".
                          format(1.0/(wnumber[i]*pc.MTC), isoID[i], elow[i]), 2)
      i += 1


    # Calculate the wavelength in microns:
    wlength[:] = 1.0 / (wnumber * pc.MTC)
    # Set isotopic index to start counting from 0:
    isoID -= 1
    isoID[np.where(isoID < 0)] = 9 # 10th isotope had index 0 --> 10-1=9
    # Calculate gf:
    Ia = np.asarray(self.isoratio)
    gf = (S0 * pc.C1 * PFzero[isoID] /  Ia[isoID] *
             np.exp( pc.C2 * elow / self.T0)  /
                 (1.0-np.exp(-pc.C2 * wnumber    / self.T0)) )
    # Alternative way:
    gf2 = A21 * g2 * pc.C1 / (8.0 * np.pi * pc.c) / wnumber**2.0

    # FINDME: Delete me when gf-vs-gf2 issue solved.
    # print(gf)
    # print(gf2)
    # print((gf/gf2)[:20])
    # print(isoID[:20])
    # print(np.amax(gf/gf2), np.amin(gf/gf2))
    # import matplotlib.pyplot as plt
    # plt.figure(1)
    # plt.plot(isoID, gf/gf2, "b,")
    # plt.ylim(0,4)
    # plt.xlim(-0.5, 5.5)
    # plt.savefig("gf.png")
    # plt.figure(2)
    # print(np.shape(self.gi), np.shape(PFzero))
    # ggi = np.asarray(self.gi)[isoID]
    # plt.plot(ggi, gf2/gf, "b,")
    # plt.plot([0,12], [0,12], "r")
    # plt.ylim(0,15)
    # plt.xlim(-0.5, 12.5)
    # plt.savefig("gf2.png")

    pt.msg(verbose-20, "GF values: {}".format(gf))
    data.close()
    return wnumber, gf2, elow, isoID
