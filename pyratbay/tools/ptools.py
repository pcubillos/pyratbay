# Copyright (c) 2015-2016 Patricio Cubillos and contributors.
# Pyrat Bay is open-source software under the FINDME license.

__all__ = ["parray", "msg", "warning", "error",
           "binsearch", "pprint", "divisors", "u", "unpack", "sep"]

import sys, os
import traceback
import textwrap
import struct
import numpy as np

from .. import constants as pc

"""
Tools for the Pyrat-Bay project.
"""

# Warning/error banner:
sep = 70*":"


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


def msg(verblevel, message, file=None, indent=0, si=None, noprint=False):
  """
  Conditional message printing to screen.

  Parameters
  ----------
  verblevel: Integer
     If positive, print the given message.
  message: String
     Message to print.
  file: File pointer
     If not None, print message to the given file pointer.
  indent: Integer
     Number of blank spaces for indentation.
  si: Integer
     Subsequent indentation.  If None, keep indent as the subsequent
     indentation.
  noprint: Boolean
     If True, do not print and return the string instead.
  """
  if verblevel <= 0:
    return

  # Output text message:
  text = ""

  # Indentation strings:
  indspace  = " "*indent
  sindspace = indspace
  if si is not None:
    sindspace = " "*si

  # Break the text down into sentences (line-breaks):
  sentences = message.splitlines()
  for s in sentences:
    line = textwrap.fill(s, break_long_words=False, break_on_hyphens=False,
                         initial_indent=indspace, subsequent_indent=sindspace)
    text += line + "\n"

  # Do not print, just return the string:
  if noprint:
    return text

  # Print to screen:
  print(text[:-1])  # Remove the trailing "\n"
  sys.stdout.flush()
  # Print to file, if requested:
  if file is not None:
    file.write(text)
    file.flush()


def warning(message, wlog, file=None):
  """
  Print message surrounded by colon banners.
  Append message to wlog.
  Add message to file if not None.

  Parameters
  ----------
  message: String
     Message to print.
  wlog:  list of strings
     List of warning messages.
  file: File pointer
     If not None, also print to the given file.
  """
  # Wrap the message:
  text = msg(1, message, indent=4, noprint=True)[:-1]
  # Add banners around:
  warntext = "\n{:s}\n  Warning:\n{:s}\n{:s}\n".format(sep, text, sep)

  # Append warning message to warnings log:
  wlog.append(text)
  # Print warning message to screen:
  print(warntext)
  sys.stdout.flush()
  # Print warning message to file (if requested):
  if file is not None:
    file.write(warntext + "\n")
    file.flush()


def error(message, file=None, lev=-2):
  """
  Pretty print error message.

  Parameters
  ----------
  message: String
     Message to print.
  file: File pointer
     If not None, also print to the given file.
  """
  # Trace back the file, function, and line where the error source:
  t = traceback.extract_stack()
  # Extract fields:
  modpath    = t[lev][0]                       # Module path
  modname    = modpath[modpath.rfind('/')+1:]  # Madule name
  funcname   = t[lev][2]                       # Function name
  linenumber = t[lev][1]
  # Text to print:
  text = ("\n{:s}\n  Error in module: '{:s}', function: '{:s}', line: {:d}\n"
          "{:s}\n{:s}".format(sep, modname, funcname, linenumber,
                              msg(1,message,indent=4,noprint=True)[:-1], sep))

  # Print to screen:
  print(text)
  sys.stdout.flush()
  # Print to file (if requested):
  if file is not None:
    file.write(text)
    file.close()
  sys.exit(0)


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
  Index of record with ...
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


def divisors(number):
  """
  Find all the integer divisors of number.
  """
  divs = []
  for i in np.arange(1, number/2+1):
    if number % i == 0:
      divs.append(i)
  divs.append(number)
  return np.asarray(divs, np.int)


def unpack(file, n, dtype):
  """
  Wrapper for struct unpack.

  Parameters:
  -----------
  file: File object
     File object to read from.
  n: Integer
     Number of elements to read from file.
  dtype: String
     Data type of the bytes read.

  Returns:
  --------
  output: Scalar, tuple, or string
     If dtype is 's' return the string.
     If there is a single element to read, return the scalar value.
     Else, return a tuple with the elements read.
  """
  # Compute the reading format:
  fmt  = "{:d}{:s}".format(n, dtype)
  # Calculate the number of bytes to read:
  size = struct.calcsize(fmt)
  # Read:
  output = struct.unpack(fmt, file.read(size))
  # Return:
  if (n == 1) or (dtype == "s"):
    return output[0]
  else:
    return output


def u(units):
  """
  Get the conversion factor (to the CGS system) for units.

  Parameters:
  -----------
  units: String
     Name of units
  """
  # Accept only valid units:
  if units not in ["eV",  "A", "nm", "um", "mm", "cm", "m", "km", "au", "pc",
                   "mbar", "bar", "kelvin", "amu", "amagat"]:
    # Throw error:
    error("Units name '{:s}' does not exist.".format(units), lev=-3)
  exec("factor = pc.{:s}".format(units))
  return factor

