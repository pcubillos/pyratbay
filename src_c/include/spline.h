// Copyright (c) 2021 Patricio Cubillos
// Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

/* Interpolate cubic spline. See Section 3.3 of Numerical Recipes. */

void cubic_spline(
        PyArrayObject *xin,   /* Input X array                          */
        PyArrayObject *yin,   /* Input Y array                          */
        PyArrayObject *y2nd,  /* Second derivatives of yi at xi         */
        PyArrayObject *xout,  /* X array for interpolated values        */
        PyArrayObject *yout,  /* Output (spline interpolated) array     */
        int lo,               /* Lower index of xout to interpolate     */
        int hi){              /* Upper index of xout to interpolate     */

    int i=0, n, nin;
    double dx, a, b;

    nin = PyArray_DIM(xin, 0);

    for (n=lo; n<=hi; n++){
        /* Do a binary search: */
        i = binsearchapprox(xin, INDd(xout,n), i, nin-1);
        /* Enforce: xin[i] <= xout[n] (except if xin[i] == xin[nin-1]): */
        if (i == nin-1 || INDd(xout,n) < INDd(xin,i))
            i--;

        dx = INDd(xin,(i+1)) - INDd(xin,i);
        a = (INDd(xin,(i+1)) - INDd(xout,n)) / dx;
        b = (INDd(xout,n)-INDd(xin,i)) / dx;

        INDd(yout,n) = a*INDd(yin,i) + b*INDd(yin,(i+1))
            + ((a*a*a-a)*INDd(y2nd,i) + (b*b*b-b)*INDd(y2nd,(i+1))) * dx*dx/6.0;
    }
    return;
}
