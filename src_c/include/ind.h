// Copyright (c) 2021-2022 Patricio Cubillos
// Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

// Definitions for indexing Numpy arrays:

// 1D double ndarray:
#define INDd(a,i) *((double *)(PyArray_DATA(a) + i * PyArray_STRIDE(a, 0)))

// 2D double ndarray:
#define IND2d(a,i,j) *((double *)( \
    PyArray_DATA(a) \
    + i * PyArray_STRIDE(a, 0) \
    + j * PyArray_STRIDE(a, 1)))

// 3D double ndarray:
#define IND3d(a,i,j,k) *((double *)( \
    PyArray_DATA(a) \
    + i * PyArray_STRIDE(a, 0) \
    + j * PyArray_STRIDE(a, 1) \
    + k * PyArray_STRIDE(a, 2)))

// 4D double ndarray:
#define IND4d(a,i,j,k,l) *((double *)( \
    PyArray_DATA(a) \
    + i * PyArray_STRIDE(a, 0) \
    + j * PyArray_STRIDE(a, 1) \
    + k * PyArray_STRIDE(a, 2) \
    + l * PyArray_STRIDE(a, 3)))

// 1D integer ndarray:
#define INDi(a,i) *((int *)(PyArray_DATA(a) + i * PyArray_STRIDE(a, 0)))

// 2D integer ndarray:
#define IND2i(a,i,j) *((int *)( \
    PyArray_DATA(a) \
    + i * PyArray_STRIDE(a, 0) \
    + j * PyArray_STRIDE(a, 1)))
