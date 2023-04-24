// Copyright (c) 2021-2023 Patricio Cubillos
// Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <math.h>

#include "ind.h"
#include "utils.h"
#include "spline.h"


PyDoc_STRVAR(second_deriv__doc__,
"Compute the second derivative of an array y(x).  \n\
                                                  \n\
Parameters                                        \n\
----------                                        \n\
yin: 1D double ndarray                            \n\
    Input y values.                               \n\
xin: 1D double ndarray                            \n\
    Input x array.                                \n\
                                                  \n\
Return                                            \n\
------                                            \n\
y2nd: 1D double ndarray                           \n\
    Second derivatives of yin at xin.");

static PyObject *second_deriv(PyObject *self, PyObject *args){
    PyArrayObject *xin, *yin, *y2nd;
    double sig, p, *u;
    npy_intp nin[1];
    int i, n;

    /* Load inputs: */
    if (!PyArg_ParseTuple(args, "OO", &yin, &xin))
        return NULL;

    /* Allocate arrays: */
    nin[0] = PyArray_DIM(yin,0);
    y2nd = (PyArrayObject *) PyArray_SimpleNew(1, nin, NPY_DOUBLE);

    n = (int)nin[0] - 1;
    u  = calloc(n, sizeof(double));

    INDd(y2nd,0) = INDd(y2nd,n) = 0.0;
    u[0] = 0.0;
    /* Calculate second derivatives: */
    for (i=1; i<n; i++){
        sig = (INDd(xin,i) - INDd(xin,(i-1))) /
              (INDd(xin,(i+1)) - INDd(yin,(i-1)));
        p = sig * INDd(y2nd,(i-1)) + 2.0;
        INDd(y2nd,i) = (sig-1.0)/p;
        u[i] = (INDd(yin,(i+1))-INDd(yin,i)) / (INDd(xin,(i+1))-INDd(xin,i))
            - (INDd(yin,i)-INDd(yin,(i-1))) / (INDd(xin,i)-INDd(xin,(i-1)));
        u[i] = (6.0*u[i] / (INDd(xin,(i+1))-INDd(xin,(i-1)))-sig*u[i-1])/p;
    }

    for (i=n-1; i>=0; i--){
        INDd(y2nd,i) = INDd(y2nd,i)*INDd(y2nd,(i+1)) + u[i];
    }

    /* Free arrays: */
    free(u);
    return Py_BuildValue("N", y2nd);
}


PyDoc_STRVAR(splinterp_1D__doc__,
"Cubic spline interpolation.                                         \n\
                                                                     \n\
Parameters                                                           \n\
----------                                                           \n\
yin: 1D double ndarray                                               \n\
    Input y values.                                                  \n\
xin: 1D double ndarray                                               \n\
    Input x array.                                                   \n\
y2nd: 2D double ndarray                                              \n\
    Second derivatives of yin at xin.                                \n\
xout: 1D double ndarray                                              \n\
    X array for inpterpolated values.                                \n\
extrap: double                                                       \n\
    Value assigned to extrapolated points.                           \n\
                                                                     \n\
Return                                                               \n\
------                                                               \n\
yout: 1D float ndarray                                               \n\
    Spline interpolated values.");


static PyObject *splinterp_1D(PyObject *self, PyObject *args){
    PyArrayObject *xin, *yin, *y2nd, *xout, *yout;
    npy_intp nout[1];
    int nin, epleft, epright;
    double extrap;

    /* Load inputs: */
    if (!PyArg_ParseTuple(args, "OOOOd", &yin, &xin, &y2nd, &xout, &extrap))
        return NULL;

    /* Get the number of datapoints: */
    nin = (int)PyArray_DIM(xin,  0);
    nout[0] = PyArray_DIM(xout, 0);

    /* Allocate arrays: */
    yout = (PyArrayObject *) PyArray_SimpleNew(1, nout, NPY_DOUBLE);

    /* Extrapolated values: */
    epleft = 0;
    while(INDd(xout, epleft)  < INDd(xin,0)){
        INDd(yout, epleft) = extrap;
        epleft += 1;
    }
    epright = (int)nout[0]-1;
    while(INDd(xout, epright) > INDd(xin,(nin-1))){
        INDd(yout, epright) = extrap;
        epright -= 1;
    }

    /* Interpolate: */
    cubic_spline(xin, yin, y2nd, xout, yout, epleft, epright);

    return Py_BuildValue("N", yout);
}


PyDoc_STRVAR(splinterp_2D__doc__,
"Cubic spline interpolation for 2D array along the first axis,       \n\
Repeating across second axis.                                        \n\
                                                                     \n\
Parameters                                                           \n\
----------                                                           \n\
yin: 2D double ndarray                                               \n\
    Input y values of shape [nin, nin2].                             \n\
xin: 1D double ndarray                                               \n\
    Input x array of length nin.                                     \n\
z: 2D double ndarray                                                 \n\
    Second derivatives of yin at xin (of shape [nin, nin2]).         \n\
xout: 1D double ndarray                                              \n\
    X value where to interpolate (of length nout).                   \n\
yout: 2D double ndarray                                              \n\
    X value where to interpolate (of shape [nout,nin2]).             \n\
lo: Integer                                                          \n\
    Lower index along second axis of yin to interpolate.             \n\
hi: Integer                                                          \n\
    Upper index along second axis of yin to interpolate.");

static PyObject *splinterp_2D(PyObject *self, PyObject *args){
    PyArrayObject *xin, *yin, *z, *xout, *yout;
    int i, j, nin, nout, index, lo, hi;
    double dx, dy, deltax, a, b, c;

    /* Load inputs: */
    if (!PyArg_ParseTuple(
            args,
            "OOOOOii",
            &yin, &xin, &z, &xout, &yout, &lo, &hi))
        return NULL;

    nin = (int)PyArray_DIM(xin, 0);
    nout = (int)PyArray_DIM(xout, 0);

    for (i=0; i<nout; i++){
        if (INDd(xout,i) < INDd(xin,0) || INDd(xout,i) > INDd(xin,(nin-1)))
            return Py_BuildValue("d", NAN);

        /* Binary search to find index: */
        index = binsearchapprox(xin, INDd(xout,i), 0, nin-1);
        /* Enforce: x[i] <= xout (except if x[N-1] == xout): */
        if (index == nin-1 || INDd(xout,i) < INDd(xin,index)){
            index--;
        }

        /* If the requested x falls on a given x, no interpolation needed: */
        if (INDd(xin,index) == INDd(xout,i)){
            for (j=lo; j<hi; j++){
                IND2d(yout,i,j) = IND2d(yin,index,j);
            }
            continue;
        }
        /* Else, spline-interpolate: */
        dx = INDd(xin,(index+1)) - INDd(xin,index);
        deltax = INDd(xout,i) - INDd(xin,index);
        for (j=lo; j<hi; j++){
            dy = IND2d(yin,(index+1),j) - IND2d(yin,index,j);
            a = (IND2d(z,(index+1),j) - IND2d(z,index,j)) / (6*dx);
            b = 0.5 * IND2d(z,index,j);
            c = dy/dx - dx/6 * (IND2d(z,(index+1),j) + 2*IND2d(z,index,j));
            IND2d(yout,i,j) = IND2d(yin,index,j)
                              + deltax*(c + deltax*(b + deltax*a));
        }
    }

    return Py_BuildValue("d", 0.0);
}


PyDoc_STRVAR(lin_interp_2D__doc__,
"Linear interpolation for 2D array along the first axis,             \n\
Repeating across second axis.                                        \n\
                                                                     \n\
Parameters                                                           \n\
----------                                                           \n\
yin: 2D double ndarray                                               \n\
    Input y values of shape [nin, nin2].                             \n\
xin: 1D double ndarray                                               \n\
    Input x array of length nin.                                     \n\
dy_dx: 2D double ndarray                                             \n\
    Derivatives of yin at xin (of shape [nin-1, nin2-1]).            \n\
xout: 1D double ndarray                                              \n\
    X value where to interpolate (of length nout).                   \n\
yout: 2D double ndarray                                              \n\
    X value where to interpolate (of shape [nout,nin2]).             \n\
lo: Integer                                                          \n\
    Lower index along second axis of yin to interpolate.             \n\
hi: Integer                                                          \n\
    Upper index along second axis of yin to interpolate.");

static PyObject *lin_interp_2D(PyObject *self, PyObject *args){
    PyArrayObject *xin, *yin, *dy_dx, *xout, *yout;
    int i, j, nin, nout, index, lo, hi;
    double deltax;

    // Load inputs:
    if (!PyArg_ParseTuple(
            args,
            "OOOOOii",
            &yin, &xin, &dy_dx, &xout, &yout, &lo, &hi))
        return NULL;

    nin = (int)PyArray_DIM(xin, 0);
    nout = (int)PyArray_DIM(xout, 0);

    for (i=0; i<nout; i++){
        if (INDd(xout,i) < INDd(xin,0) || INDd(xout,i) > INDd(xin,(nin-1)))
            return Py_BuildValue("d", NAN);

        // Binary search to find index:
        index = binsearchapprox(xin, INDd(xout,i), 0, nin-1);
        // Enforce: x[i] <= xout (except if x[N-1] == xout):
        if (index == nin-1 || INDd(xout,i) < INDd(xin,index)){
            index--;
        }

        // If the requested x falls on a given x, no interpolation needed:
        if (INDd(xin,index) == INDd(xout,i)){
            for (j=lo; j<hi; j++){
                IND2d(yout,i,j) = IND2d(yin,index,j);
            }
            continue;
        }
        // Else, spline-interpolate
        deltax = INDd(xout,i) - INDd(xin,index);
        for (j=lo; j<hi; j++){
            IND2d(yout,i,j) = IND2d(yin,index,j) + deltax*IND2d(dy_dx,index,j);
        }
    }

    return Py_BuildValue("d", 0.0);
}


/* The module doc string */
PyDoc_STRVAR(
    spline__doc__,
    "Python wrapper for Cubic Spline Interpolation."
);


/* A list of all the methods defined by this module. */
static PyMethodDef spline_methods[] = {
    {"second_deriv", second_deriv, METH_VARARGS, second_deriv__doc__},
    {"splinterp_1D", splinterp_1D, METH_VARARGS, splinterp_1D__doc__},
    {"splinterp_2D", splinterp_2D, METH_VARARGS, splinterp_2D__doc__},
    {"lin_interp_2D", lin_interp_2D, METH_VARARGS, lin_interp_2D__doc__},
    {NULL, NULL, 0, NULL}
};


/* Module definition for Python 3. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_spline",
    spline__doc__,
    -1,
    spline_methods
};

/* When Python 3 imports a C module named 'X' it loads the module */
/* then looks for a method named "PyInit_"+X and calls it. */
PyObject *PyInit__spline (void) {
    PyObject *module = PyModule_Create(&moduledef);
    import_array();
    return module;
}
