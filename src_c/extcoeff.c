// Copyright (c) 2016 Patricio Cubillos and contributors.
// Pyrat Bay is currently proprietary software (see LICENSE).

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <stdarg.h>

#include "ind.h"
#include "constants.h"
#include "utils.h"

PyDoc_STRVAR(extinction__doc__,
"Calculate the extinction-coefficient for each molecule, at\n\
a given pressure and temperature, over a wavenumber range. \n\
                                                           \n\
Parameters                                                 \n\
----------                                                 \n\
ext: 2D float ndarray                                      \n\
   Output extinction coefficient where to put the results. \n\
profile: 1D float ndarray                                  \n\
   Array of Voigt profiles.                                \n\
psize: 2D integer ndarray                                  \n\
   Profiles half-size.                                     \n\
pindex: 2D integer ndarray                                 \n\
   Index where each profile starts.                        \n\
lorentz: 1D Float ndarray                                  \n\
   Sample of Lorentz HWHMs.                                \n\
doppler: 1D Float ndarray                                  \n\
   Sample of Doppler HWHMs.                                \n\
wn: 1D Float ndarray                                       \n\
   Spectrum wavenumber array (cm-1).                       \n\
own: 1D Float ndarray                                      \n\
   Oversampled wavenumber array (cm-1).                    \n\
divisors: 1D integer ndarray                               \n\
   Integer divisors for oversampling factor.               \n\
moldensity: 1D Float ndarray                               \n\
   Atmospheric species number density (molecules cm-3).    \n\
molq: 1D Float ndarray                                     \n\
   Atmospheric species mole mixing ratio.                  \n\
molrad: 1D Float ndarray                                   \n\
   Atmospheric species collision radius (A).               \n\
molmass: 1D Float ndarray                                  \n\
   Atmospheric species mass (gr mol-1).                    \n\
isoimol: 1D Float ndarray                                  \n\
   Isotopes species index (from atmospheric-species array).\n\
isomass: 1D Float ndarray                                  \n\
   Isotopes mass (gr mol-1).                               \n\
isoratio: 1D Float ndarray                                 \n\
   Isotopes abundance ratio.                               \n\
isoz: 1D Float ndarray                                     \n\
   Isotopes partition function.                            \n\
isoiext: 1D Float ndarray                                  \n\
   Isotopes species index in ext. coeff. table.            \n\
lwn: 1D Float ndarray                                      \n\
   Line-transition wavenumber (cm-1).                      \n\
elow: 1D Float ndarray                                     \n\
   Line-transition lower-state energy (cm-1).              \n\
gf: 1D Float ndarray                                       \n\
   Line-transition oscillator strength.                    \n\
lID: 1D integer ndarray                                    \n\
   Line-transition isotope ID.                             \n\
ethresh: Float                                             \n\
   Extinction-coefficient threshold factor.                \n\
pressure: Float                                            \n\
   Atmospheric-layer pressure (barye).                     \n\
temp: Float                                                \n\
   Atmospheric-layer temperature (K).                      \n\
logtext: String                                            \n\
   Char array where to store log text outputs.             \n\
verb: Integer                                              \n\
   Verbosity level.                                        \n\
add: Integer                                               \n\
   Flag, if True calculate the extinction coefficient (in cm-1) \n\
   for this layer, if False calculate the extinction coefficient\n\
   (in cm2 molecules-1) for each species in the layer.          \n\
                                                                \n\
Uncredited developers                                           \n\
---------------------                                           \n\
Patricio Rojo  U. Cornell.");

static PyObject *extinction(PyObject *self, PyObject *args){
  PyArrayObject *ext,                             /* Extinction coefficient */
           *profile, *psize, *pindex,
           *lorentz, *doppler,                    /* Voigt data             */
           *wn, *own, *divisors,                  /* Wavenumber data        */
           *moldensity, *molq, *molrad, *molmass, /* Molecular data         */
           *isoimol, *isomass, *isoratio,
           *isoz, *isoiext,                       /* Isotopic data          */
           *lwn, *lID, *elow, *gf;                /* Line transition data   */

  long nmol, niso, nlines, next, ndivs,
       onwn, dnwn,
       minj, maxj,
       nLor, nDop;
  int iwn, iown, idwn, offset, subw,
      imol, ofactor, iprof,
      nadd=0, nskip=0, neval=0,
      verb, add=0;
  int i, j, m, ln; /* Auxilliary for-loop indices    */
  long jj;
  double pressure, temp, csdiameter, density, minwidth=1e5, vwidth, ethresh,
         florentz, fdoppler, wnstep, ownstep, dwnstep, wavn, next_wn, k;
  double *alphal, *alphad, *kmax, **ktmp, *kprop;
  int *idop, *ilor;
  char *logtext;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOdddsi|i",
                              &ext,
                              &profile, &psize, &pindex, &lorentz, &doppler,
                              &wn, &own, &divisors,
                              &moldensity, &molq, &molrad, &molmass,
                              &isoimol, &isomass, &isoratio, &isoz, &isoiext,
                              &lwn, &elow, &gf, &lID,
                              &ethresh, &pressure, &temp,
                              &logtext, &verb, &add))
    return NULL;

  nLor   = PyArray_DIM(lorentz,  0);  /* Number of Lorentz widths            */
  nDop   = PyArray_DIM(doppler,  0);  /* Number of Doppler widths            */
  nmol   = PyArray_DIM(molmass,  0);  /* Number of species                   */
  niso   = PyArray_DIM(isomass,  0);  /* Number of isotopes                  */
  ndivs  = PyArray_DIM(divisors, 0);  /* Number of divisors of osamp         */
  onwn   = PyArray_DIM(own,      0);  /* Number of fine-wavenumber samples   */
  nlines = PyArray_DIM(lwn,      0);  /* Number of line transitions          */
  next   = PyArray_DIM(ext,      0);  /* Number of extinction-coef. species  */

  if (add)
    next = 1;

  strcpy(logtext, "\0");
  /* Constant factors for line widths:                                      */
  fdoppler = sqrt(2*KB*temp/AMU) * SQRTLN2 / LS;
  florentz = sqrt(2*KB*temp/PI/AMU) / LS;

  /* Allocate alpha Lorentz and Doppler arrays:                             */
  alphal = (double *)malloc(niso * sizeof(double));
  alphad = (double *)malloc(niso * sizeof(double));

  /* Allocate width indices array:                                          */
  idop = (int *)malloc(niso * sizeof(int));
  ilor = (int *)malloc(niso * sizeof(int));

  /* Allocate line strength per transition:                                 */
  kprop = (double *)malloc(nlines * sizeof(double));
  /* Array to hold maximum line strength:                                   */
  kmax  = (double *)calloc(next,    sizeof(double));


  ktmp    = (double **)malloc(next *     sizeof(double *));
  ktmp[0] = (double  *)calloc(next*onwn, sizeof(double  ));
  for (i=1; i<next; i++)
    ktmp[i] = ktmp[0] + onwn*i;

  /* Calculate the isotopes' widths:                                        */
  for (i=0; i<niso; i++){
    imol = INDi(isoimol, i);
    /* Lorentz profile width:                                               */
    alphal[i] = 0.0;
    for(j=0; j<nmol; j++){
      /* Isotope's collision diameter with j-th species:                    */
      csdiameter = INDd(molrad, imol) + INDd(molrad, j);
      /* Number density of molecule colliding with current isotope:         */
      density = INDd(molq,j) * pressure / (KB*temp);
      /* Line width:                                                        */
      alphal[i] += density * csdiameter * csdiameter *
                   sqrt(1/INDd(isomass,i) + 1/INDd(molmass,j));
      if (i==-1)
        printf("j:%d,  %.6e,  %.6e,  %.6e\n", j,
                     density, csdiameter, INDd(molmass,j));
    }
    alphal[i] *= florentz;

    /* Doppler profile width (divided by central wavenumber):               */
    alphad[i] = fdoppler / sqrt(INDd(isomass,i));
    /* Print Lorentz and Doppler broadening widths:                         */
    if(i <= 0){
      //printf("Imass: %.4f, imol: %d\n", INDd(isomass,i), imol);
      //printf("Len %d\n", n);
      msg(verb-6, logtext, "    Lorentz: %.3e cm-1, Doppler: %.3e cm-1.\n",
          alphal[i], alphad[i]*INDd(wn,0), temp, pressure);
      //printf("Len %d\n", n);
    }
    /* Estimate the Voigt width:                                            */
    vwidth = 0.5346*alphal[i] + sqrt(pow(alphal[i], 2)*0.2166    +
                                     pow(alphad[i]*INDd(wn,0),2) );
    minwidth = fmin(minwidth, vwidth);

    /* Search for aDop and aLor indices for alphal[i] and alphad[i]:        */
    idop[i] = binsearchapprox(doppler, alphad[i]*INDd(wn,0), 0, (int)nDop-1);
    ilor[i] = binsearchapprox(lorentz, alphal[i],            0, (int)nLor-1);
  }

  wnstep  = INDd(wn, 1) - INDd(wn, 0);
  ownstep = INDd(own,1) - INDd(own,0);
  /* Set the wavenumber sampling resolution:                                */
  for (i=1; i<ndivs; i++)
    if (INDi(divisors,i)*ownstep >= 0.5 * minwidth){
      break;
    }
  ofactor = INDi(divisors,(i-1));   /* Dynamic-sampling oversampling factor */
  dwnstep = ownstep * ofactor;      /* Dynamic-sampling grid stepsize       */
  dnwn    = 1 + (onwn-1) / ofactor; /* Number of dynamic-sampling values    */
  msg(verb-6, logtext, "    Dynamic-sampling grid interval: %.4e "
              "(factor:%i, minwidth:%.3e)\n", dwnstep, ofactor, minwidth);
  msg(verb-6, logtext, "    Number of dynamic-sampling values: %ld\n", dnwn);

  /* Find the maximum line-strength per molecule:                           */
  for (ln=0; ln<nlines; ln++){
    wavn = INDd(lwn, ln);     /* Wavenumber of line transition              */
    i    = INDi(lID, ln);     /* Isotope index of line transition:          */
    m    = INDi(isoiext, i);  /* Extinction-coefficient table index         */
    if (add)
      m = 0;  /* Collapse everything into the first index                   */

    /* If this line falls beyond limits, skip to next line transition:      */
    if ((wavn < INDd(own,0)) || (wavn > INDd(own,(onwn-1))))
      continue;

    /* Calculate the line strength divided by the molecular abundance:      */
    kprop[ln] = k = INDd(isoratio,i)       *  /* Density fraction           */
        SIGCTE * INDd(gf,ln)               *  /* Constant * gf              */
        exp(-EXPCTE*INDd(elow,ln) / temp)  *  /* Level population           */
        (1-exp(-EXPCTE*wavn/temp))         /  /* Induced emission           */
        INDd(isoz,i);                         /* Partition function         */
    /* Check if this is the maximum k:                                      */
    kmax[m] = fmax(kmax[m], k);
  }

  /* Compute the extinction-coefficient for each species:                   */
  for (ln=0; ln<nlines; ln++){
    wavn = INDd(lwn, ln);    /* Wavenumber                                  */
    i    = INDi(lID, ln);    /* Isotope index                               */
    if (add)
      m  = 0;                /* First index of extinction-coefficient       */
    else
      m  = INDi(isoiext, i); /* Extinction-coefficient table index          */

    if ((wavn < INDd(own,0)) || (wavn > INDd(own,(onwn-1))))
      continue;

    /* Line strength:                                                       */
    k = kprop[ln];

    /* Index of closest oversampled wavenumber:                             */
    iown = (wavn - INDd(own,0))/ownstep;
    if (fabs(wavn - INDd(own,(iown+1))) < fabs(wavn - INDd(own,iown)))
      iown++;

    /* Check if the next line falls on the same sampling index:             */
    while (ln+1 != nlines && INDi(lID, (ln+1)) == i){
      next_wn = INDd(lwn, (ln+1));
      if (fabs(next_wn - INDd(own,iown)) < ownstep){
        nadd++;
        ln++;
        /* Add the contribution from this line into the opacity:            */
        k += kprop[ln];
      }
      else
        break;
    }

    /* Estimate the line-transition width:                                  */
    vwidth = 0.5346*alphal[i] + sqrt(pow(alphal[i], 2.0)*0.2166 +
                                     pow(alphad[i]*wavn, 2.0));

    /* Index of wavenumber closest and smaller than wavn:                   */
    iwn = (wavn - INDd(wn,0))/wnstep;
    /* Line is more than 10 Voigt HW from nearest wn sampling division:     */
    if (fabs(wavn-INDd(wn,iwn)-0.5*wnstep) > 10.0*vwidth){
      if (wavn - INDd(wn,iwn) > 0.5*wnstep)
        iwn++;

      if (add)  /* Multiply by the species' density:                        */
        k *= INDd(moldensity, (INDi(isoimol,i)));
      IND2d(ext, m, iwn) += k/wnstep;
      continue;
    }

    /* Skip weakly contributing lines:                                      */
    if (k < ethresh * kmax[m]){
      nskip++;
      continue;
    }

    /* Multiply by the species' density:                                    */
    if (add)
      k *= INDd(moldensity, (INDi(isoimol,i)));

    /* Index of closest (but not larger than) dynamic-sampling wavenumber:  */
    idwn = (wavn - INDd(own,0))/dwnstep;

    /* Calculate index for Doppler width:                                   */
    //if (alphad[i]*wavn/alphal[i] >= 0.1)
    idop[i] = pyramidsearch(doppler, alphad[i]*wavn, idop[i], (int)nDop-1);

    /* Sub-sampling offset between center of line and dyn-sampled wn:       */
    subw = iown - idwn*ofactor;
    /* Offset between the profile and the wavenumber-array indices:         */
    offset = ofactor*idwn - IND2i(psize, idop[i], ilor[i]) + subw;
    /* Range that contributes to the opacity:                               */
    /* Set the lower and upper indices of the profile to be used:           */
    minj = idwn - (IND2i(psize, idop[i], ilor[i]) - subw) / ofactor;
    maxj = idwn + (IND2i(psize, idop[i], ilor[i]) + subw) / ofactor;
    if (minj < 0)
      minj = 0;
    if (maxj > dnwn)
      maxj = dnwn;

    /* Add the contribution from this line to the opacity spectrum:         */
    iprof = IND2i(pindex, idop[i], ilor[i]);
    jj = iprof + ofactor*minj - offset;
    for(j=(int)minj; j<maxj; j++){
      ktmp[m][j] += k * INDd(profile, jj);
      jj += ofactor;
    }
    neval++;
  }
  msg(verb-5, logtext, "    Number of co-added lines:     %8i  (%5.2f%%)\n",
              nadd,  nadd*100.0/nlines);
  msg(verb-5, logtext, "    Number of skipped profiles:   %8i  (%5.2f%%)\n",
              nskip, nskip*100.0/nlines);
  msg(verb-5, logtext, "    Number of evaluated profiles: %8i  (%5.2f%%)\n",
              neval, neval*100.0/nlines);

  /* Downsample ktmp to the final sampling size:                            */
  for (m=0; m<next; m++){
    downsample(ktmp, ext, (int)dnwn, (int)round(wnstep/ownstep/ofactor), m);
  }

  /* Free the no-longer used memory:                                        */
  free(alphal);
  free(alphad);
  free(kmax);
  free(idop);
  free(ilor);
  free(ktmp[0]);
  free(ktmp);
  free(kprop);

  return Py_BuildValue("i", 1);
}


PyDoc_STRVAR(interp_ec__doc__,
"Interpolate the extinction coefficient.                            \n\
                                                                    \n\
Parameters:                                                         \n\
-----------                                                         \n\
extinction: 3D float ndarray                                        \n\
   Extinction coefficient array [nwave] to calculate.               \n\
etable 1D float ndarray                                             \n\
   Tabulated extinction coefficient [nmol, ntemp, nwave].           \n\
ttable: 1D float ndarray                                            \n\
   Tabulated temperature array [ntemp].                             \n\
mtable: 1D float ndarray                                            \n\
   Tabulated molID of species [nmol].                               \n\
temperature: Float                                                  \n\
   Atmospheric layer temperature.                                   \n\
density: 1D float ndarray                                           \n\
   Density of species in the atmospheric layer [nmol].              \n\
molID: 1D integer ndarray                                           \n\
   Molecular ID of species in the atmosphere [nmol].");

static PyObject *interp_ec(PyObject *self, PyObject *args){
  PyArrayObject *extinction,  /* Interpolated extinction coefficient        */
                *etable,      /* Tabulated extinction coefficient           */
                *ttable,      /* Tabulated temperature                      */
                *mtable,      /* Tabulated species ID                       */
                *density,     /* Density of species at given layer          */
                *molID;       /* mol ID of the atmospheric species          */

  int tlo, thi, imol;
  long nwave, nmol, nmolID, ntemp;  /* Number of wavenumber, species, temps */
  int i, j;                /* Auxilliary for-loop indices                   */
  double ext,              /* Temporary extinction coefficient              */
         temperature;      /* Atmospheric-layer temperature                 */

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOOdOO",  &extinction, &etable,
                                          &ttable, &mtable,
                                          &temperature,
                                          &density, &molID))
    return NULL;
  nmol   = PyArray_DIM(etable, 0);  /* Number of species samples            */
  ntemp  = PyArray_DIM(etable, 1);  /* Number of temperature samples        */
  nwave  = PyArray_DIM(etable, 2);  /* Number of spectral samples           */
  nmolID = PyArray_DIM(molID,  0);  /* Number of species in atmosphere      */

  /* Find index of grid-temperature immediately lower than temperature:     */
  tlo = binsearchapprox(ttable, temperature, 0, (int)ntemp-1);
  if (temperature < INDd(ttable,tlo)){
    tlo--;
  }
  thi = tlo + 1;

  /* Add contribution from each molecule:                                   */
  for (j=0; j<nmol; j++){
    imol = valueinarray(molID, INDi(mtable,j), (int)nmolID);
    for (i=0;  i<nwave; i++){
      /* Linear interpolation of the extinction coefficient:                */
      ext = (IND3d(etable,j,tlo,i) * (INDd(ttable,thi) - temperature) +
             IND3d(etable,j,thi,i) * (temperature - INDd(ttable,tlo)) ) /
            (INDd(ttable,thi) - INDd(ttable,tlo));

      INDd(extinction, i) += INDd(density,imol) * ext;
    }
  }

  return Py_BuildValue("i", 1);
}


/* The module doc string    */
PyDoc_STRVAR(extcoeff__doc__, "Python wrapper for the extinction-coefficient "
                              "calculation.");

/* A list of all the methods defined by this module.                        */
static PyMethodDef extcoeff_methods[] = {
    {"extinction", extinction, METH_VARARGS, extinction__doc__},
    {"interp_ec",  interp_ec,  METH_VARARGS, interp_ec__doc__},
    {NULL,         NULL,       0,            NULL}              /* sentinel */
};


#if PY_MAJOR_VERSION >= 3
/* Module definition for Python 3.                                          */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "extcoeff",
    extcoeff__doc__,
    -1,
    extcoeff_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit_extcoeff (void) {
  PyObject *module = PyModule_Create(&moduledef);
  import_array();
  return module;
}

#else
/* When Python 2 imports a C module named 'X' it loads the module           */
/* then looks for a method named "init"+X and calls it.                     */
void initextcoeff(void){
  Py_InitModule3("extcoeff", extcoeff_methods, extcoeff__doc__);
  import_array();
}
#endif
