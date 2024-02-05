// Copyright (c) 2021-2024 Patricio Cubillos
// Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "constants.h"
#include "utils.h"


PyDoc_STRVAR(
    alkali_cross_section__doc__,
"TBD                                                             \n\
                                                                 \n\
Parameters                                                       \n\
----------                                                       \n\
wn: 1D float ndarray                                             \n\
    Wavenumber spectrum (cm-1).                                  \n\
temp: 1D float ndarray                                           \n\
    Temperature at each layer (K).                               \n\
TBD                                                              \n\
                                                                 \n\
Returns                                                          \n\
-------                                                          \n\
TBD                                                              \n\
");

static PyObject *alkali_cross_section(PyObject *self, PyObject *args){
    PyArrayObject *wn, *temp, *pressure, *wn0, *gf, *dwave, *voigt_det,
        *i_wn0, *ec;
    int i, j, k, t, flip, nwave, nlayers, nlines;
    double detuning_wn, dsigma, mass, lorentz_par, part_func, cutoff,
        lorentz, doppler, full_width, dwn, abs_dwn;

    /* Load inputs: */
    if (!PyArg_ParseTuple(
            args,
            "OOOOOdddddOOOO",
            &pressure, &wn, &temp, &voigt_det, &ec,
            &detuning_wn, &mass, &lorentz_par, &part_func, &cutoff,
            &wn0, &gf, &dwave, &i_wn0
    ))
        return NULL;

    /* Get the spectrum size: */
    nlayers = (int)PyArray_DIM(temp, 0);
    nwave = (int)PyArray_DIM(wn, 0);
    nlines = (int)PyArray_DIM(wn0, 0);

    // if wn is in decreasing order, flip the indexing
    flip = signbit(INDd(wn,1)-INDd(wn,0));

    // Calculate the cross sections:
    for (j=0; j<nlines; j++){
        for (i=0; i<nlayers; i++){
            // Doppler half width (cm-1):
            doppler =
                sqrt(2.0 * KB * INDd(temp,i) / (mass*AMU))
                * INDd(wn0,j) / LS;
            // Lorentz half width (cm-1):
            lorentz =
                lorentz_par
                * pow(INDd(temp,i)/2000.0,-0.7) * INDd(pressure,i) / ATM;

            full_width = 2.0 * (
                0.5346*lorentz +
                sqrt(0.2166*pow(lorentz,2.0) + pow(doppler,2.0))
            );
            // Detuning frequency (cm-1):
            dsigma = detuning_wn * pow(INDd(temp,i)/500.0, 0.6);

            for (k=0; k<nwave; k++){
                if (flip)
                    t = nwave - k - 1;
                else
                    t = k;
                dwn = INDd(wn,t) - INDd(wn0,j);
                abs_dwn = fabs(dwn);
                if (dwn < -cutoff)
                    continue;
                else if (dwn > cutoff)
                    break;

                // Extinction outside the detuning region (power law):
                if (abs_dwn >= dsigma)
                    IND2d(ec,i,t) +=
                        IND2d(voigt_det,i,j) * pow(abs_dwn/dsigma,-1.5)
                        * C3 * INDd(gf,j) / part_func
                        * exp(-C2*(abs_dwn-dsigma) / INDd(temp,i));
                else {
                    // Lorentz profile in the core (detuning region):
                    IND2d(ec,i,t) +=
                        lorentz / PI / (pow(lorentz,2.0) + pow(dwn,2.0)) *
                        C3 * INDd(gf,j) / part_func;
                    // Note this equation neglects exp(-Elow/T)*(1-exp(-wn0/T))
                    // because it is approximately 1.0 at T=[100K--4000K]
                }
            }
        }
    }
    return Py_BuildValue("i", 1);
}



// The module doc string
PyDoc_STRVAR(
    alkali__doc__,
    "Wrapper for the alkali extinction-coeficient calculation."
);

// A list of all the methods defined by this module:
static PyMethodDef alkali_methods[] = {
    {"alkali_cross_section", alkali_cross_section, METH_VARARGS, alkali_cross_section__doc__},
    {NULL, NULL, 0, NULL}  // sentinel
};


// Module definition for Python 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_alkali",
    alkali__doc__,
    -1,
    alkali_methods
};

// When Python 3 imports a C module named 'X' it loads the module
// then looks for a method named "PyInit_"+X and calls it.
PyObject *PyInit__alkali (void) {
    PyObject *module = PyModule_Create(&moduledef);
    import_array();
    return module;
}
