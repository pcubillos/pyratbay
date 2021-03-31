// Copyright (c) 2021 Patricio Cubillos
// Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)


/* Do the final sum for the Simpson integration algorithm.  Based on
   the Python implementation:
   github.com/scipy/scipy/blob/v0.15.1/scipy/integrate/quadrature.py */
double
simpson(PyArrayObject *y,         /* Values of function to integrate        */
        PyArrayObject *hsum,      /* See simpson.c's geth function          */
        PyArrayObject *hratio,    /* See simpson.c's geth function          */
        PyArrayObject *hfactor,   /* See simpson.c's geth function          */
        int n){                   /* Number of elements in y                */

    int i=0, j;
    double res = 0.0;
    /* Add contribution from each interval: */
    for (i=0; i < (n-1)/2; i++){
        j = 2*i;
        res += (INDd(y, (j  )) * (2.0 - 1.0/INDd(hratio, i)) +
                INDd(y, (j+1)) * INDd(hfactor, i)            +
                INDd(y, (j+2)) * (2.0 - INDd(hratio, i)) ) * INDd(hsum,i);
    }
    return res/6.0;
}


/* Same as above, but for a 2D input function */
double
simpson2(PyArrayObject *y,
        PyArrayObject *hsum,
        PyArrayObject *hratio,
        PyArrayObject *hfactor,
        int n,
        int k){  /* Index at which to integrate */

    int i=0, j;
    double res=0.0;
    /* Add contribution from each interval: */
    for (i=0; i < (n-1)/2; i++){
        j = 2*i;
        res += (IND2d(y,  j,    k) * (2 - 1.0/INDd(hratio, i)) +
                IND2d(y, (j+1), k) * INDd(hfactor, i)          +
                IND2d(y, (j+2), k) * (2 - INDd(hratio, i)) ) * INDd(hsum,i);
    }
    return res/6.0;
}

