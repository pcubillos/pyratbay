// Copyright (c) 2016-2019 Patricio Cubillos and contributors.
// Pyrat Bay is currently proprietary software (see LICENSE).

double
simpson(PyArrayObject *y,         /* Values of function to integrate        */
        PyArrayObject *hsum,      /* See simpson.c's geth function          */
        PyArrayObject *hratio,    /* See simpson.c's geth function          */
        PyArrayObject *hfactor,   /* See simpson.c's geth function          */
        int n){                   /* Number of elements in y                */

  /* Do the final sum for the Simpson integration algorithm.  Based on 
     the Python implementation:
     github.com/scipy/scipy/blob/v0.15.1/scipy/integrate/quadrature.py      */

  int i=0,           /* for-loop index                                      */
      j;             /* Array index for each interval                       */
  double res = 0.0;  /* The results                                         */
  int even=1-(n%2);
  /* Add contribution from each interval:                                   */
  for (i=0; i < (n-1)/2; i++){
    /* Skip first value of y if there's an even number of samples:          */
    j = 2*i + even;
    res += (INDd(y, (j  )) * (2.0 - INDd(hratio, i))     +
            INDd(y, (j+1)) * INDd(hfactor, i)            +
            INDd(y, (j+2)) * (2.0 - 1.0/INDd(hratio, i)) ) * INDd(hsum,i);
  }

  return res/6.0;
}
