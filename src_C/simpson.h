double
simpson(PyArrayObject *y,      PyArrayObject *hsum,
        PyArrayObject *hratio, PyArrayObject *hfactor, int n){
  /* Do the final sum for the Simpson integration algorithm.  Based on 
     the Python implementation:
     github.com/scipy/scipy/blob/v0.15.1/scipy/integrate/quadrature.py      */
  int i=0, j;
  double res = 0.0;
  /* Last case, value limited between consecutive indices of array:         */
  for (i=0; i<n; i++){
    /* Skip first interval if there are even samples:                       */
    j = 2*i + (n%2==1);
    res += (INDd(y, (j  )) * (2.0 - INDd(hratio, i))     +
            INDd(y, (j+1)) * INDd(hfactor, i)            +
            INDd(y, (j+2)) * (2.0 - 1.0/INDd(hratio, i)) ) * INDd(hsum,i);
  }

  return res/6.0;
}
