/* Solve the Natural-Cubic-Spline Tridiagonal System.
   See: Kincaid & Cheney (2008): Numerical Mathematics & Computing,
        Algorithm 1 (page 391).                                             */
void
tri(PyArrayObject *z,    /* Second derivatives of y                         */
    PyArrayObject *y,    /* y values of array to interpolate from           */
    PyArrayObject *dx,   /* Distance between consecutive x values           */
    int n){              /* Number of elements in y                         */

  int i;  /* Auxiliary for-loop indices                                     */
  double *b, *u, *v;
  b = calloc(n-1, sizeof(double));
  u = calloc(n-1, sizeof(double));
  v = calloc(n-1, sizeof(double));

  /* Step 1:                                                                */
  for (i=0; i<n-1; i++){
    b[i] = (INDd(y,(i+1)) - INDd(y,i)) / INDd(dx,i);
  }

  /* Step 2:                                                                */
  u[1] = 2 * (INDd(dx,1) + INDd(dx,0));
  v[1] = 6 * (b[1]       - b[0]);

  for (i=2; i<n-1; i++){
    u[i] = 2*(INDd(dx,i)+INDd(dx,(i-1))) - INDd(dx,(i-1))*INDd(dx,(i-1))/u[i-1];
    v[i] = 6*(b[i]      - b[i-1]       ) - v[i-1]        *INDd(dx,(i-1))/u[i-1];
  }

  /* Step 3:                                                                */
  INDd(z,0) = INDd(z,(n-1)) = 0;  /* First and last are zeroes              */
  for (i=n-2; i>0; i--){
    INDd(z,i) = (v[i] - INDd(dx,i) * INDd(z,(i+1))) / u[i];
  }

  free(b);
  free(u);
  free(v);
  return;
}


/* Evaluate the cubic spline  y(x) interpolating from yi(xi).
   See: Kincaid & Cheney (2008): Numerical Mathematics & Computing,
        Equation 12, (page 392).                                            */
void
spline3(PyArrayObject *xin,   /* Input X array                              */
        PyArrayObject *yin,   /* Input Y array                              */
        int nin,              /* Length of xin                              */
        PyArrayObject *xout,  /* X array for interpolated values            */
        PyArrayObject *yout,  /* Output (spline interpolated) array         */
        int ninterp,          /* Number of points to interpolate            */
        PyArrayObject *z,     /* Second derivatives of yi at xi             */
        PyArrayObject *dx,    /* Spacing between xin points                 */
        int nep){             /* Number of extrapolated points              */

  int i=0, n;  /* Indices                                                   */
  double B;

  /* Calculate the spline interpolation y-value for each x-coordinate
     in the requested array                                                 */
  for (n=nep; n<nep+ninterp; n++){
    /* Else, do a binary search:                                            */
    i = binsearchapprox(xin, INDd(xout,n), i, nin-1);
    /* Enforce: xin[i] <= xout[n] (except if xin[i] == xin[nin-1]):         */
    if (i == nin-1 || INDd(xout,n) < INDd(xin,i)){
      i--;
    }

    /* Factor for linear coefficient:                                       */
    B = (INDd(yin,(i+1)) - INDd(yin,i)) / INDd(dx,i) -
        INDd(dx,i)/6 * (INDd(z,(i+1)) + 2 * INDd(z,i));

    /* Calculate y-value from cubic polynomial:                             */
    INDd(yout,n) = INDd(yin,i) + (INDd(xout,n) - INDd(xin,i)) * B +
      pow(INDd(xout,n)-INDd(xin,i),2)*0.5*INDd(z,i)  +
      pow(INDd(xout,n)-INDd(xin,i),3)*(INDd(z,(i+1))-INDd(z,i))/(6*INDd(dx,i));
  }
  return;
}
