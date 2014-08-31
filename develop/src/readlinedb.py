import struct, sys
import numpy as np
import objects as o
import ptools  as pt
import pconstants as pc

def readlinedb(pyrat):
  """
  Read line transition data.

  Modification History:
  ---------------------
  2014-06-08  patricio  Implemented readheader.
  2014-06-15  patricio  Implemented checkrange.
  2014-08-15  patricio  Implemented setimol.
  """
  pt.msg(pyrat.verb, "\nReading line transition info:")

  # Count number of TLI files:
  pyrat.lt.nTLI = len(pyrat.linedb)
  TLI = []
  TLIiline = [] # unused

  # Read data bases header info:
  for n in np.arange(pyrat.lt.nTLI):
    # Open-read TLI data base:
    TLI.append(open(pyrat.linedb[n], "r"))
    # Read headers info:
    TLIiline.append(readheader(pyrat, TLI[n]))

  # Set link to molecule's indices:
  setimol(pyrat) 

  # Read line-transition data:
  for n in np.arange(pyrat.lt.nTLI):
    readlinetransition(pyrat, TLI[n], n)

  pt.msg(pyrat.verb, "Done.\n")


def readheader(pyrat, linefile):
  """
  Open and read database headers info.

  Returns:
  --------
  iline: Integer
  Pointer position in linefile where the line-transition info begins.

  Modification History:
  ---------------------
  2014-06-08  patricio  Initial implementation.
  2014-08-03  patricio  Updated to pylineread5.0.
  2014-08-17  patricio  Store molname in DB object.
  """
  # Read header:
  endian = linefile.read(1)
  if endian == 'b':
    pt.msg(pyrat.verb, "TLI data storage: Big endian",    2)
  elif endian == 'l':
    pt.msg(pyrat.verb, "TLI data storage: Little endian", 2)

  # Compare endianness:
  endianness = sys.byteorder
  if endianness[0:1] != endian:
    pt.error("Incompatible endianness between TLI file (%s) and pyrat "
             "(%s)."%(endian, endianness[0:1]))

  # Read TLI version:
  TLI_ver, LR_ver, LR_rev = struct.unpack("3h", linefile.read(6))
  pt.msg(pyrat.verb, "TLI version: %s  Lineread version: %s  "
                     "Lineread revision: %s"%(TLI_ver, LR_ver, LR_rev), 2)

  # Read initial and final wavenumber from TLI:
  lt_wni, lt_wnf = struct.unpack("2d", linefile.read(16))
  pt.msg(pyrat.verb, "TLI wavenumber range (cm^-1): [%.1f, %.1f]"%(
                                                       lt_wni, lt_wnf), 2)
  # Check TLI and pyrat wavelength ranges:
  checkrange(pyrat, lt_wni, lt_wnf)

  # Read number of data bases: 
  Ndb = struct.unpack("h", linefile.read(2))[0]
  pt.msg(pyrat.verb, "Number of data bases: %d"%Ndb, 2)

  # Cumulative isotope index:
  acumiso = pyrat.iso.niso

  for i in np.arange(Ndb):
    db = o.database()
    # Read Database name:
    lenDBname = struct.unpack("h",             linefile.read(2))[0]
    db.name   = struct.unpack("%ds"%lenDBname, linefile.read(lenDBname))[0]
    pt.msg(pyrat.verb, "Data base name: '{:s}'".format(db.name), 2)
    # Read Molecule name:
    lenMolec  = struct.unpack("h",             linefile.read(2))[0]
    db.molname   = struct.unpack("%ds"%lenMolec, linefile.read(lenMolec))[0]
    pt.msg(pyrat.verb, "Molecule name: '{:s}'".format(db.molname), 2)
    # Read temperature array:
    db.ntemp, db.niso = struct.unpack("hh", linefile.read(4))
    db.temp = np.asarray(struct.unpack("%dd"%db.ntemp,
                                       linefile.read(8*db.ntemp)))
    pt.msg(pyrat.verb, 'Temperature range: {:4.0f} -- {:4.0f} K.'.format(
                        db.temp[0], db.temp[-1]), 2)

    # Allocate arrays for isotopic info:
    name    = np.zeros(db.niso, pc.strfmt)
    mass    = np.zeros(db.niso)
    ratio   = np.zeros(db.niso)
    dbindex = np.zeros(db.niso, np.int)
    db.z    = np.zeros((db.niso, db.ntemp))
    #db.c    = np.zeros((db.niso, db.ntemp))

    # Store per-isotope info:    
    for j in np.arange(db.niso):
      dbindex[j] = i
      lenIsoName = struct.unpack("h",              linefile.read(2))[0]
      name[j]    = struct.unpack("%ds"%lenIsoName, linefile.read(lenIsoName))[0]
      mass[j]    = struct.unpack("d",              linefile.read(8))[0]
      ratio[j]   = struct.unpack("d",              linefile.read(8))[0]
      db.z[j]    = np.asarray(struct.unpack("%dd"%db.ntemp,
                                            linefile.read(8*db.ntemp)))
      #db.c[j]    = np.asarray(struct.unpack("%dd"%db.ntemp,
      #                                      linefile.read(8*db.ntemp)))
      # Give info:
      pt.msg(pyrat.verb, "Isotope: %s,  mass: %.4f,  ratio: %.5g"%(
                         name[j], mass[j], ratio[j]), 3)
      pt.msg(pyrat.verb, "Z = [%.2e, %.2e, ..., %.2e]"%(db.z[j,0],
                                                 db.z[j, 1], db.z[j,-1]), 4)
      #pt.msg(pyrat.verb, "C = [%.2e, %.2e, ..., %.2e]"%(db.c[j,0],
      #                                           db.c[j, 1], db.c[j,-1]), 3)

    # Add the number of isotopes read:
    pyrat.iso.niso += db.niso
    pt.msg(pyrat.verb-10, "pyrat niso: {:d}".format(pyrat.iso.niso), 2)

    # Store name, mass in isotopes structure:
    pyrat.iso.name    = np.concatenate((pyrat.iso.name,    name))
    pyrat.iso.mass    = np.concatenate((pyrat.iso.mass,    mass))
    pyrat.iso.dbindex = np.concatenate((pyrat.iso.dbindex, dbindex))

    # Set isotope correlative index for DB:
    db.iiso  = acumiso
    acumiso += db.niso
    pyrat.lt.db.append(db)
    #cdb, cui = struct.unpack("hh", linefile.read(4))
    pt.msg(pyrat.verb, "DB index: {:d}, Cumulative Iso: {:d}".format(i,
                                                                   acumiso), 2)

  # Keep count of number of databases:
  pyrat.lt.ndb  += Ndb

  # Return the pointer position in lineinfo file:
  return linefile.tell()


def readlinetransition(pyrat, linefile, dbindex):
  """
  Read the databases line transition info.

  Modification History:
  ---------------------
  2014-06-22  patricio  Initial implementation
  """
  # Read the number of line transitions:
  nTransitions = struct.unpack("l", linefile.read(8))[0]
  # Position where the line-transition data begins:
  initrec = linefile.tell()
  # Count the number of transitions:
  linefile.seek(0, 2)
  endrec = linefile.tell()
  nrec = (endrec - initrec + 0.0) / pc.tlireclen
  if nrec != nTransitions:
    pt.error("The remaining data file size (%.1f) does not correspond "
             "to the number of transitions (%d)."%(nrec, nTransitions))
  pt.msg(pyrat.verb, "There are %d line transitions in TLI "
                     "file."%nTransitions, 2)

  # Show boundary values:
  linefile.seek(initrec)
  ini = struct.unpack('d', linefile.read(8))[0]
  linefile.seek((nrec-1)*pc.dreclen + initrec, 0)
  fin = struct.unpack('d', linefile.read(8))[0]
  pt.msg(pyrat.verb, "Database wavenumber boundaries: %.3f -- %.3f cm^-1"%(
                     ini, fin), 2)

  # pyrat low wavelength boundary in microns:
  #pwl_low  = 1.0/(pyrat.wnhigh*pc.units['um'])
  #pt.msg(pyrat.verb, "Targer initial wavelength: %.7f um"%pwl_low, 2)
  #pyrat.wnlow = 3448.0
  start = pt.binsearch(linefile, pyrat.wnlow, initrec, nTransitions-1, False)
  linefile.seek(start*pc.dreclen + initrec, 0)
  pt.msg(pyrat.verb, "Found initial transition (%d):  %.7f cm^-1"%(start,
                     struct.unpack('d', linefile.read(8))[0]), 2)

  #pwl_high = 1.0/(pyrat.wnlow*pc.units['um'])
  #pt.msg(pyrat.verb, "Targer final wavelength:   %.7f um"%pwl_high, 2)
  #pyrat.wnhigh = 3571.48
  end = pt.binsearch(linefile, pyrat.wnhigh, initrec, nTransitions-1, True)
  linefile.seek(end*pc.dreclen + initrec, 0)
  pt.msg(pyrat.verb, "Found final transition (%d):    %.7f cm^-1"%(end,
                     struct.unpack('d', linefile.read(8))[0]), 2)

  # Number of transitions to read:
  nread = end - start + 1
  pt.msg(pyrat.verb, "Reading %d transitions."%nread, 2)

  # Allocate arrays:
  wn    = np.zeros(nread) 
  gf    = np.zeros(nread) 
  elow  = np.zeros(nread) 
  isoid = np.zeros(nread, np.short)

  # Read data into arrays:
  linefile.seek(start*pc.dreclen + initrec, 0)
  wn[:]    = struct.unpack('%dd'%nread, linefile.read(nread*pc.dreclen))
  linefile.seek(start*pc.sreclen + nTransitions*pc.dreclen + initrec, 0)
  isoid[:] = struct.unpack('%dh'%nread, linefile.read(nread*pc.sreclen))
  linefile.seek(start*pc.dreclen + nTransitions*(pc.dreclen+pc.sreclen) +
                initrec, 0)
  elow[:]  = struct.unpack('%dd'%nread, linefile.read(nread*pc.dreclen))
  linefile.seek(start*pc.dreclen + nTransitions*(2*pc.dreclen+pc.sreclen) +
                initrec, 0)
  gf[:]    = struct.unpack('%dd'%nread, linefile.read(nread*pc.dreclen))

  # for i in np.arange(nread):
  #   wn[i], isoid[i] = struct.unpack('dh', linefile.read(10))
  #   elow[i], gf[i]  = struct.unpack('dd', linefile.read(16))

  # FINDME: Need number of isotopes per TLI file (or Number of DBs)
  #isoid += 0*pyrat.lt.db[dbindex].iiso  # Add isotope correlative index
  # FINDME: Need to sort transition  if there is more than 1 TLI file

  pyrat.lt.wn    = np.concatenate((pyrat.lt.wn,    wn)   )
  pyrat.lt.gf    = np.concatenate((pyrat.lt.gf,    gf)   )
  pyrat.lt.elow  = np.concatenate((pyrat.lt.elow,  elow) )
  pyrat.lt.isoid = np.concatenate((pyrat.lt.isoid, isoid))
  pyrat.lt.ntransitions += nread

  #print(pyrat.lt.db[dbindex].iiso)


def checkrange(pyrat, wn_low, wn_high):
  """
  Display a warning if line database spectral range does not completely
  include the pyrat spectral range.

  Parameters:
  -----------
  pyrat: Object
  wn_low: Float
    Database's wavenumber lower boundary (in cm^-1).
  wn_high: Float
    Database's wavenumber higher boundary (in cm^-1).

  Modification History:
  ---------------------
  2014-06-15  patricio  Initial implementation
  2015-08-03  patricio  Changed wavelength to wavenumber
  """
  # Print out warning if ranges dont overlap:
  if (wn_low > pyrat.wnhigh) or (wn_high < pyrat.wnlow):
    pt.warning("TLI wavenumber range (%.2f - %.2f cm^-1) doesn't overlap with "
               "Pyrat wavenumber range (%.2f - %.2f cm^-1)."%(wn_low, wn_high,
                pyrat.wnlow, pyrat.wnhigh))
  # Print out warning if TLI range is smaller than the pyrat required range:
  elif (wn_low > pyrat.wnlow) or  (wn_high < pyrat.wnhigh):
    pt.warning("TLI wavenumber range (%.2f - %.2f cm^-1) doesn't contain the "
               "full Pyrat wavenumber range (%.2f - %.2f cm^-1)."%(wn_low,
               wn_high, pyrat.wnlow, pyrat.wnhigh))


def setimol(pyrat):
  """
  Set the molecule index for the list of isotopes.

  Modification History:
  ---------------------
  2014-08-17  patricio  Initial implemetation.
  """
  # Allocate imol array:
  pyrat.iso.imol = np.zeros(pyrat.iso.niso, np.int)
  # For each isotope:
  for i in np.arange(pyrat.iso.niso):
    # Get molecule name from database object:
    molname = pyrat.lt.db[pyrat.iso.dbindex[i]].molname
    # Set index:
    if molname not in pyrat.mol.symbol:  # Isotope's not in molecule list
      pyrat.iso.imol[i] = -1
    else:
      pyrat.iso.imol[i] = np.where(pyrat.mol.symbol == molname)[0]
  pt.msg(pyrat.verb, "Isotope's molecule indices:\n"
                     "  {:s}".format(str(pyrat.iso.imol)), 2)

