#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "expn.h"

double
xi(double gamma, double tau){
  return 2.0/3.0 * ( (1/gamma) * (1 + (0.5*gamma*tau-1) * exp(-gamma*tau))   +
                     gamma * (1 - 0.5*pow(tau,2)) * expn(2, gamma*tau) + 1.0 );
}


PyDoc_STRVAR(ptfunc__doc__,
"Generats a PT profile based on input free parameters and pressure   \n\
array.  If no inputs are provided, it will run in demo mode, using   \n\
free parameters given by the Line 2013 paper and some dummy pressure \n\
parameters.                                                          \n\
                                                                     \n\
Implementation from equations (13)-(16) of Line et al. (2013), Apj,  \n\
  775, 137.                                                          \n\
                                                                     \n\
Inputs                                                               \n\
------                                                               \n\
params: 1D float ndarray                                             \n\
  Array of free parameters:                                          \n\
   log10(kappa):  Planck thermal IR opacity in units cm^2/gr         \n\
   log10(gamma1): Visible-to-thermal stream Planck mean opacity ratio\n\
   log10(gamma2): Visible-to-thermal stream Planck mean opacity ratio\n\
   alpha:         Visible-stream partition (0.0--1.0)                \n\
   beta:          A 'catch-all' for albedo, emissivity, and day-night\n\
                  redistribution (on the order of unity)             \n\
pressure: 1D float ndarray                                           \n\
  Array of pressure values in bars.                                  \n\
temperature: 1D float ndarray                                        \n\
  Output temperature array in Kelvin.                                \n\
R_star: Float                                                        \n\
   Stellar radius (in meters).                                       \n\
T_star: Float                                                        \n\
   Stellar effective temperature (in Kelvin degrees).                \n\
T_int:  Float                                                        \n\
   Planetary internal heat flux (in Kelvin degrees).                 \n\
sma:    Float                                                        \n\
   Semi-major axis (in meters).                                      \n\
grav:   Float                                                        \n\
   Planetary surface gravity in cm/second^2.                         \n\
                                                                     \n\
  Returns                                                            \n\
  -------                                                            \n\
  T: temperature array                                               \n\
                                                                     \n\
  Example:                                                           \n\
  --------                                                           \n\
  >>> import PT as pt                                                \n\
  >>> import scipy.constants as sc                                   \n\
  >>> import matplotlib.pyplot as plt                                \n\
  >>> import numpy as np                                             \n\
                                                                     \n\
  >>> Rsun = 6.995e8 # Sun radius in meters                          \n\
                                                                     \n\
  >>> # Pressure array (bars):                                       \n\
  >>> p = np.logspace(2, -5, 100)                                    \n\
                                                                     \n\
  >>> # Physical (fixed for each planet) parameters:                 \n\
  >>> Ts = 5040.0        # K                                         \n\
  >>> Ti =  100.0        # K                                         \n\
  >>> a  = 0.031 * sc.au # m                                         \n\
  >>> Rs = 0.756 * Rsun  # m                                         \n\
  >>> g  = 2192.8        # cm s-2                                    \n\
                                                                     \n\
  >>> # Fitting parameters:                                          \n\
  >>> kappa  = -1.5   # log10(3e-2)                                  \n\
  >>> gamma1 = -0.8   # log10(0.158)                                 \n\
  >>> gamma2 = -0.8   # log10(0.158)                                 \n\
  >>> alpha  = 0.5                                                   \n\
  >>> beta   = 1.0                                                   \n\
  >>> params = [kappa, gamma1, gamma2, alpha, beta]                  \n\
  >>> T0 = pt.PT(p, params, Rs, Ts, Ti, a, g)                        \n\
                                                                     \n\
  >>> plt.figure(1)                                                  \n\
  >>> plt.clf()                                                      \n\
  >>> plt.semilogy(T0, p, lw=2, color='b')                           \n\
  >>> plt.ylim(p[0], p[-1])                                          \n\
  >>> plt.xlim(800, 2000)                                            \n\
  >>> plt.xlabel('Temperature  (K)')                                 \n\
  >>> plt.ylabel('Pressure  (bars)')                                 \n\
                                                                     \n\
  Developers:                                                        \n\
  -----------                                                        \n\
  Madison Stemm      astromaddie@gmail.com                           \n\
  Patricio Cubillos  pcubillos@fulbrightmail.org");

static PyObject *pt(PyObject *self, PyObject *args){
  PyArrayObject *freepars, *pressure, *temperature;
  double Rstar, Tstar, Tint, sma, grav,
         kappa, gamma1, gamma2, alpha, beta, tau, Tirr, xi1, xi2;
  int i, nlayers;     /* Auxilliary for-loop indices                        */

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOddddd", &freepars, &pressure, &temperature,
                                          &Rstar, &Tstar, &Tint, &sma, &grav))
    return NULL;

  /* Get array size:                                                        */
  nlayers = PyArray_DIM(pressure, 0);

  /* Unpack the model parameters:                                           */
  kappa  = pow(10.0, INDd(freepars, 0));
  gamma1 = pow(10.0, INDd(freepars, 1));
  gamma2 = pow(10.0, INDd(freepars, 2));
  alpha  = INDd(freepars, 3);
  beta   = INDd(freepars, 4);

  /* Stellar input temperature (at top of atmosphere):                      */
  Tirr = beta * sqrt(Rstar / (2.0*sma)) * Tstar;

  /* Gray IR optical depth:                                                 */
  for (i=0; i<nlayers; i++){
    tau = kappa * (INDd(pressure,i)*1e6) / grav;  /* bars to barye          */
    xi1 = xi(gamma1, tau);
    xi2 = xi(gamma2, tau);

    INDd(temperature,i) = pow(0.75 * (pow(Tint,4) * (2.0/3.0 + tau) +
                                      pow(Tirr,4) * (1-alpha) * xi1 +
                                      pow(Tirr,4) *    alpha  * xi2 ), 0.25);
  }

  return Py_BuildValue("O", temperature);
}

/* The module doc string                                                    */
PyDoc_STRVAR(pt__doc__,
  "Python wrapper for Total Internal Partition Sum (TIPS) calculation.");

/* A list of all the methods defined by this module.                        */
static PyMethodDef pt_methods[] = {
    {"pt",      pt,      METH_VARARGS, ptfunc__doc__},
    {NULL,      NULL,    0,            NULL}
};


/* When Python imports a C module named 'X' it loads the module */
/* then looks for a method named "init"+X and calls it.         */
void initpt(void){
  Py_InitModule3("pt", pt_methods, pt__doc__);
  import_array();
}
