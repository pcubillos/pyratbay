import sys, os, traceback, textwrap, struct
import numpy as np
#from mpi4py import MPI

import pconstants as pc
"""
pyrat tools: Tools for pyrat project.

"""

def parray(string):
  """
  Convert a string containin a list of white-space-separated (and/or
  newline-separated) values into a numpy array
  """
  if string == 'None':
    return None
  try:    # If they can be converted into doubles, do it:
    return np.asarray(string.split(), np.double)
  except: # Else, return a string array:
    return string.split()


def exit(comm=None, abort=False, message=None, comm2=None):
  """
  Stop execution.

  Parameters:
  -----------
  comm: MPI communicator
     An MPI Intracommunicator.
  abort: Boolean
     If True send (gather) an abort flag integer through comm.
  message: String
     Print message on exit.

  Modification History:
  ---------------------
  2014-04-20  patricio  Initial implementation (extracted from transit project).
  """
  if message is not None:
    print(message)
  if comm is not None:
    if abort:
      #comm_gather(comm, np.array([1], dtype='i'), MPI.INT)
      pass
    comm.Barrier()
    comm.Disconnect()
  if comm2 is not None:
    comm2.Barrier()
    comm2.Disconnect()
  sys.exit(0)


def error(message, lev=-2):
  """
  Pretty print error message

  2015-01-18  patricio  Added lev argument to print.
  """
  # Trace back the file, function, and line where the error source:
  t = traceback.extract_stack()
  # Extract fields:
  efile = t[lev][0]
  efile = efile[efile.rfind('/')+1:]
  efunc = t[lev][2]
  eline = t[lev][1]
  # Indent and wrap message to 70 characters:
  msg = textwrap.fill(message, initial_indent   ="    ",
                               subsequent_indent="    ")
  # Print it out:
  print(
    "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
    "  Error in module: '%s', function: '%s', line: %d\n"
    "%s\n"
    "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"%
    (efile, efunc, eline, msg))
  sys.exit(0)


def warning(message):
  """
  Print message surrounded by colon bands.

  Modification History:
  ---------------------
  2014-06-15  patricio  Initial implementation.
  """
  print(
    "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
    "  Warning:")
  msg(1, message, 4)
  print(
    "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n")


def msg(verblevel, message, indent=0):
  """
  Conditional message printing to screen.

  Modification History:
  ---------------------
  2014-06-15  patricio  Added Documentation.
  """
  sentences = message.splitlines()
  indspace = " "*indent
  if verblevel > 0:
    for s in sentences:
      msg = textwrap.fill(s, replace_whitespace=True,
                          initial_indent=indspace, subsequent_indent=indspace)
      print(msg)


def binsearch(dbfile, wavelength, rec0, nrec, upper=True):
  """
  Do a binary search in TLI dbfile for record that has wavelength iwl

  Parameters:
  -----------
  dbfile: File object
     TLI file where to search.
  wavelength: Scalar
     Target wavelength.
  rec0: Integer
     Position of first record.
  nrec: Integer
     Number of records.
  upper: Boolean
     Wavelength is an upper boundary, return the index of value smaller
     than wavelength.

  Returns:
  --------
  Index of record with

  Modification History:
  ---------------------
  2013        madison   Initial implementation.
  2014-03-05  patricio  Added documentation, updated Madison's code, and
                        included searchup parameter. pcubillos@fulbrightmail.org
  2014-06-22  patricio  Adapted from pylineread project to read from TLI.
  2014-08-03  patricio  Added boundaries condition in sequential search.
  """
  # Wavelength of record:
  ilo = 0
  ihi = nrec

  # Start binary search:
  while ihi - ilo > 1:
    # Middle record index:
    irec = (ihi + ilo)/2

    # Read wavelength from TLI file record:
    dbfile.seek(irec*pc.dreclen + rec0, 0)
    rec_wl = struct.unpack('d', dbfile.read(8))[0]
    # Update search limits:
    if rec_wl > wavelength:
      ihi = irec
    else:
      ilo = irec
  #print("Binary found:      %.7f (%d)"%(rec_wl, irec))

  # Sequential search:
  if rec_wl < wavelength:
    jump =  1 # Move up
  else:
    jump = -1 # Move down

  dbfile.seek(irec*pc.dreclen + rec0, 0)
  next_wl = struct.unpack('d', dbfile.read(8))[0]
  while (np.sign(next_wl-wavelength) == np.sign(rec_wl-wavelength)):
    irec += jump
    # Check for boundaries:
    if (irec < 0 or irec > nrec):
      return np.clip(irec, 0, nrec)
    dbfile.seek(irec*pc.dreclen + rec0, 0)
    rec_wl  = next_wl
    next_wl = struct.unpack('d', dbfile.read(8))[0]
  #print("Sequential found:  %.7f"%next_wl)

  # Return the index withing the boundaries:
  return irec - (jump+upper)/2


def pprint(array, precision=3, fmt=None):
  """
  Pretty print a Numpy array.  Set desired precision and format, and
  remove line break from the string output.

  Parameters:
  -----------
  array: 1D ndarray
     Array to be pretty printed.
  precision: Integer
     Precision for floating point values.
  fmt: format
     Numeric format.

  Modification History:
  ---------------------
  2015-01-25  patricio  Initial implementation.
  """
  default_prec = np.get_printoptions().get('precision')
  np.set_printoptions(precision=precision)

  # Pretty array is a copy of array:
  parray = np.copy(array)
  if format is not None:
    parray = np.asarray(array, fmt)

  # Convert to string and remove line-breaks:
  sarray = str(array).replace("\n", "")
  np.set_printoptions(precision=default_prec)
  return sarray

