# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np

__all__ = ["eqq", "anames", "amodels"]


'''
Written by Jasmina Blecic August 3rd, 2018.
'''

def eqq(pyrat, temp, pres):
  """
  Calculate analytically species equilibrium molecular mixing ratios.
  """
  for i in np.arange(pyrat.aequil.nmodels):

    # Set the correct number of species depending on the model
    if pyrat.aequil.model[i].name == "aequil5":
      pyrat.aequil.aspecs  = ['H2O', 'CH4', 'CO', 'CO2', 'C2H2' ]

      # Calculate abundances
      pyrat.aequil.model[i].aequil5(temp, pres)

      # Allocate the abundances array
      nlayers = len(pres)
      pyrat.aequil.eqq = np.zeros((len(pyrat.aequil.aspecs), nlayers))
      pyrat.aequil.eqq += pyrat.aequil.model[i].eqq


class Aequil5():
  """
  Analytical 5 species equilibrium calculation based on Heng and Tsai (2016)
  """
  def __init__(self):
    self.name    = "aequil5"       # Model name
    self.pars    = [0.55,          # C/O        : Carbon to Oxygen elemental ratio [C/O]
                  -3.57 ]          # log10(C/H) : Carbon elemental mixing ratio 
    self.npars   = len(self.pars)  # Number of model fitting parameters
    self.parname = [r"$C/O$",      # Fitting-parameter names
                    r"$C/H$"]
    self.eqq     = None            # Equilibrium mixing fractions [Xi/sum(Xi)]

  # function to compute first equilibrium constant (K')
  def kprime(self, my_temperature, pbar):
    runiv = 8.3144621   # J/K/mol
    temperatures = np.arange(500.0, 3100.0, 100.0)
    dg = [96378.0,   72408.0,   47937.0,   23114.0,   -1949.0,  -27177.0,
         -52514.0,  -77918.0, -103361.0, -128821.0, -154282.0, -179733.0,
        -205166.0, -230576.0, -255957.0, -281308.0, -306626.0, -331911.0,
        -357162.0, -382380.0, -407564.0, -432713.0, -457830.0, -482916.0,
        -507970.0, -532995.0]
    my_dg = np.interp(my_temperature, temperatures, dg)
    result = np.exp(-my_dg/runiv/my_temperature)/pbar/pbar
    return result

  # function to compute second equilibrium constant (K2')
  def kprime2(self, my_temperature):
    runiv = 8.3144621   # J/K/mol
    temperatures = np.arange(500.0, 3100.0, 100.0)
    dg2 = [20474.0,  16689.0,  13068.0,   9593.0,   6249.0,   3021.0,   -107.0,
           -3146.0,  -6106.0,  -8998.0, -11828.0, -14600.0, -17323.0, -20000.0,
          -22634.0, -25229.0, -27789.0, -30315.0, -32809.0, -35275.0, -37712.0,
          -40123.0, -42509.0, -44872.0, -47211.0, -49528.0]
    my_dg = np.interp(my_temperature, temperatures, dg2)
    result = np.exp(-my_dg/runiv/my_temperature)
    return result

  # function to compute second equilibrium constant (K3')
  def kprime3(self, my_temperature, pbar):
    runiv = 8.3144621   # J/K/mol
    temperatures = np.arange(500.0, 3100.0, 100.0)
    dg3 = [262934.0,  237509.0,  211383.0,  184764.0,  157809.0,  130623.0, 
           103282.0,   75840.0,   48336.0,   20797.0,   -6758.0,  -34315.0, 
           -61865.0,  -89403.0, -116921.0, -144422.0, -171898.0, -199353.0,
          -226786.0, -254196.0, -281586.0, -308953.0, -336302.0, -363633.0,
          -390945.0, -418243.0]
    my_dg = np.interp(my_temperature, temperatures, dg3)
    result = np.exp(-my_dg/runiv/my_temperature)/pbar/pbar
    return result
    
  # function to compute mixing ratio for methane
  # (note: n_o is oxygen abundance, n_c is carbon abundance, kk is K')
  def n_methane(self, temp, pbar):

    self.c_o     = self.pars[0]
    self.n_c     = 10**self.pars[1]
    self.n_o     = self.n_c/self.c_o

    k1 = self.kprime(temp,pbar)
    k2 = self.kprime2(temp)
    k3 = self.kprime3(temp,pbar)
    a0 = 8.0*k1*k3*k3/k2
    a1 = 8.0*k1*k3/k2
    a2 = 2.0*k1/k2*( 1.0 + 8.0*k3*(self.n_o-self.n_c) ) + 2.0*k1*k3
    a3 = 8.0*k1/k2*(self.n_o-self.n_c) + 2.0*k3 + k1
    a4 = 8.0*k1/k2*(self.n_o-self.n_c)*(self.n_o-self.n_c) + 1.0 + 2.0*k1*(self.n_o-self.n_c)
    a5 = -2.0*self.n_c
    result = np.polynomial.polynomial.polyroots([a5, a4, a3, a2, a1, a0])
    # picks the correct root of the cubic equation
    return result[4]   

  # function to compute mixing ratio for methane
  def n_water(self, temp, pbar):
    k3     = self.kprime3(temp, pbar)
    n_ch4  = self.n_methane(temp, pbar)
    result = 2.0*k3*n_ch4*n_ch4 + n_ch4 + 2.0*(self.n_o-self.n_c)
    return result

  # function to compute mixing ratio for carbon monoxide
  def n_cmono(self, temp, pbar):
    kk     = self.kprime(temp, pbar)
    n_ch4  = self.n_methane(temp, pbar)
    n_h2o  = self.n_water(temp, pbar)
    result = kk*n_ch4*n_h2o
    return result

  # function to compute mixing ratio for carbon dioxide
  def n_cdio(self, temp, pbar):
    kk2    = self.kprime2(temp)
    n_h2o  = self.n_water(temp, pbar)
    n_co   = self.n_cmono(temp, pbar)
    result = n_co*n_h2o/kk2
    return result

  # function to compute mixing ratio for acetylene
  def n_acet(self, temp, pbar):
    kk3    = self.kprime3(temp, pbar)
    n_ch4  = self.n_methane(temp, pbar)
    result = kk3*n_ch4**2
    return result

  # calculate abundances as mixing ratio of number densities or mole numbers
  def aequil5(self, T, P):

    n_ch4  = np.zeros(len(T), dtype=float)
    n_h2o  = np.zeros(len(T), dtype=float)
    n_co   = np.zeros(len(T), dtype=float)
    n_co2  = np.zeros(len(T), dtype=float)
    n_c2h2 = np.zeros(len(T), dtype=float)
    for i in np.arange(len(T)):
      n_h2o[i]  = self.n_water(T[i], P[i])
      n_ch4[i]  = self.n_methane(T[i], P[i])
      n_co[i]   = self.n_cmono(T[i], P[i])
      n_co2[i]  = self.n_cdio(T[i], P[i])
      n_c2h2[i] = self.n_acet(T[i], P[i])

    q5 = np.vstack((n_h2o, n_ch4, n_co, n_co2, n_c2h2))
    self.eqq = q5


# List of available haze models:
amodels = [Aequil5()]

# Compile list of analytical equilibrium model names:
anames = []
for amodel in amodels:
  anames.append(amodel.name)
anames = np.asarray(anames)



