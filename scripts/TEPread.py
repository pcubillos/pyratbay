#! /usr/bin/env python

# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys
import numpy as np

class File:
  """
  Load a Transiting Extrasolar Planet (TEP) data file as a Python object,
  allowing to querry the parameters and values.

  Notes:
  ------
  The parameter 'tepfile' must be of base class FILE with READ privledges!
  This can be obtained by using the built in function 'open'.

  Input TEP files must be in the correct format:
  - The '#' character defines comments.
  - Each entry contains five fields (even if there's not a value) in this
    specific order:
      parameter name, value, uncert, units, origin/reference
    For example:
    '''
    # name       value       uncert    units                           
    planetname   HD209458b   -1        -       -                       
    startype     G0V         -1        -       SIMBAD                  
    Ts           6075        33        K       Schulder2011arXiv1103.07
    '''

  Developers:
  -----------
  Christopher Campo  UCF  ccampo@gmail.com
  Patricio Cubillos  UCF  pcubillos@fulbrightmail.org
  """
  def __init__(self, file):
    # List for the parameter and values
    self.params = np.array([])       
    self.values = []

    # Read the file
    file = open(file, 'r')
    lines = file.readlines()
    file.close()

    for line in lines:
      # Strip out comments, if any:
      try:
        line = line[0:line.index('#')].strip()
      except:
        line = line.strip()

      # Get values:
      if len(line) > 0:
        self.params = np.append(self.params, line.split()[0])
        self.values.append(np.array(line.split()[1:]))


  def evaluate(self, value):
    '''
    determines if the value is a numeric expression,
    if not, return a string.
    if it is, returns the numeric value.
    '''
    # Evaluate value if possible:
    try:
      return eval(value)
    except:
      return value


  def checkpar(self, par):
    ''' 
    Check if the input parameter exists. 
    If it does, return the reference to the values. If not, return NaN.
    '''
    try:
      id = np.where(self.params == par)[0]
      value = self.values[id]
      return value[0] if value.size == 1 else value
    except:
      return np.nan


  def getvalue(self, par, uncert=False, units=False, ref=False):
    """
    Get the values of a parameter, if it has more than one value, 
    returns an array, if not, returns the value.

    Parameters:
    -----------
    par:
    uncert: Bool
    units: Bool
    ref: Bool

    Return:
    -------
    val: If all uncert, units, and ref are False, return the parameter
         value.  Else, return a list containing the value, and the
         requested arguments.

    Examples:
    ---------
    >>> import TEPread as tr
    >>>
    >>> tep = tr.File('HD209458b.tep')
    >>> tep.getvalue('Ts')
    >>> 6075
    >>> tep.getvalue('Ts', uncert=True, units=True, ref=True)
    >>> [6075, 33, 'K', 'Schulder2011arXiv1103.0757v1']

    """
    # Returned variable:
    ret = []

    val = self.checkpar(par)
    if val is np.nan:
      print("Parameter '{:s}' not found in TEP file.".format(par))
      return np.nan

    if np.size(val) != 4:
      print("Invalid entry for '{:s}' in TEP file, there must be"
            "four fields after the parameter name.".format(par))
    # Add the parameter value:
    ret.append(self.evaluate(val[0]))
    if uncert:
      ret.append(self.evaluate(val[1]))
    if units:
      ret.append(self.evaluate(val[2]))
    if ref:
      ret.append(self.evaluate(val[3]))

    # Return list:
    if np.any([uncert, units, ref]):
      return ret
    # Return single value:
    return ret[0]


  def getstr(self, par):
    """
    Get the values of a parameter as strings.
    """
    val = self.checkpar(par)
    return val


  def pyrat_vals(self):
    """
    Print out the parameters used in the Pyrat-Bay project.
    """
    params = ["Ts", "Rs", "a", "Rp", "Mp", "loggstar"]
    print("Parameters used in the Pyrat-Bay project:")
    for i in np.arange(len(params)):
      val, unit = self.getvalue(params[i], units=True)
      print("{:8s}: {:10}  # {:s}".format(params[i], val, unit))


if __name__ == "__main__":
  """
  Open the given TEP file and print to screen the Pyrat-Bay values.
  """
  tep = File(sys.argv[1])
  tep.pyrat_vals()
