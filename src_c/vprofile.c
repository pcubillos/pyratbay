// Copyright (c) 2021-2023 Patricio Cubillos
// Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <stdarg.h>

#include "ind.h"
#include "voigt.h"
#include "utils.h"


PyDoc_STRVAR(grid__doc__,
"Calculate a set of Voigt profiles for the grid of Lorentz and Doppler     \n\
widths.                                                                    \n\
                                                                           \n\
Parameters                                                                 \n\
----------                                                                 \n\
profile: 1D double ndarray                                                 \n\
   Array (nwave) where to put the calculated Voigt profiles.               \n\
   The profiles are stacked one next to the other in the profile array.    \n\
psize: 2D integer ndarray                                                  \n\
   Array (nLor, nDop) with the half-size of spectral points of the         \n\
   profiles.  Profiles not to be calculated have psize == 0.               \n\
index: 2D integer ndarray                                                  \n\
   Array (nLor, nDop) with the index where each of the profiles start.     \n\
lorentz: 1D double ndarray                                                 \n\
   Array of Lorentz widths in cm-1.                                        \n\
doppler: 1D double ndarray                                                 \n\
   Array of Doppler widths in cm-1.                                        \n\
dwn: Float                                                                 \n\
   Wavenumber step size in cm-1.                                           \n\
verb: Integer                                                              \n\
   Verbosity flag to print info to screen.                                 \n\
                                                                           \n\
Uncredited Developers                                                      \n\
---------------------                                                      \n\
Patricio Rojo      U Cornell  pato@astro.cornell.edu (pato@oan.cl)");

static PyObject *grid(PyObject *self, PyObject *args){
    PyArrayObject *profile, *psize, *index, *doppler, *lorentz;
    double *vprofile; // Voigt profile for each (Dop,Lor) width
    double dwn;       // Wavenumber sample step size
    int nDop, nLor,  // Number of Lorentz and Doppler width samples
        nwave,   // Number of wavenumber samples of Voigt profile
        idx=0,   // Profile index position
        verb,    // Verbosity flag
        n, m, j, status;

    // Load inputs
    if (!PyArg_ParseTuple(
            args,
            "OOOOOdi",
            &profile, &psize, &index, &lorentz, &doppler, &dwn, &verb))
        return NULL;

    // Get array sizes
    nLor = (int)PyArray_DIM(lorentz, 0);
    nDop = (int)PyArray_DIM(doppler, 0);

    for (m=0; m<nLor; m++){
        for (n=0; n<nDop; n++){
            // If the profile size is > 0, calculate it
            if (IND2i(psize, m, n) != 0){
                // Number of spectral samples
                nwave = 2*IND2i(psize, m, n) + 1;
                vprofile = (double *)calloc(nwave, sizeof(double));

                if (verb>6)
                    printf(
                        "Calculating profile[%d, %d] = %d\n",
                        m, n, IND2i(psize,m,n)
                    );
                // Calculate Voigt using a width that gives an integer number
                // of 'dwn' spaced bins
                status = voigtn(
                    nwave,
                    dwn*(long)(nwave/2),
                    INDd(lorentz,m),
                    INDd(doppler,n),
                    vprofile,
                    -1,
                    nwave > _voigt_maxelements?VOIGT_QUICK:0
                );
                if (status != 1){
                    printf("voigtn() returned error code %i.\n", status);
                    return 0;
                }
                // Store values in python-object profile
                for (j=0; j<nwave; j++){
                    INDd(profile, (idx+j)) = vprofile[j];
                }
                free(vprofile);

                // Update index of profile
                IND2i(index, m, n) = idx;
                idx += 2*IND2i(psize, m, n) + 1;
            }
            else{
                // Refer to previous profile
                IND2i(index, m, n) = IND2i(index, m, (n-1));
                IND2i(psize, m, n) = IND2i(psize, m, (n-1));
                if (verb > 6)
                    printf("Skip profile[%d, %d] calculation.\n", m, n);
            }
        }
        if (verb >5)
            printf("  Calculated Voigt profile %3d/%d.\n", m+1, nLor);
    }

    return Py_BuildValue("i", 1);
}


// The module doc string    */
PyDoc_STRVAR(vprofile__doc__, "Python wrapper for Voigt profile calculation.");

// A list of all the methods defined by this module.                        */
static PyMethodDef vprofile_methods[] = {
    {"grid",   grid,   METH_VARARGS, grid__doc__},
    {NULL,     NULL,   0,            NULL}    // sentinel */
};


// Module definition for Python 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "vprofile",
    vprofile__doc__,
    -1,
    vprofile_methods
};

// When Python 3 imports a C module named 'X' it loads the module
// then looks for a method named "PyInit_"+X and calls it
PyObject *PyInit_vprofile (void) {
    PyObject *module = PyModule_Create(&moduledef);
    import_array();
    return module;
}

