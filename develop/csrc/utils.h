int
binsearchapprox(PyArrayObject *array, double value, int lo, int hi){
  /* Last case, value limited between consecutive indices of array:         */
  if (hi-lo == 1){
    /* Return closest array index to value:                                 */
    if (abs(INDd(array,hi)-value) < abs(INDd(array,lo)-value))
      return hi;
    return lo;
  }
  /* Compare to middle point and search in corresponding sub-array:         */
  else if (INDd(array,((hi+lo)/2)) > value)
    return binsearchapprox(array, value, lo, (hi+lo)/2);
  else
    return binsearchapprox(array, value, (hi+lo)/2, hi);
}


/* Downsample an array by an integer factor into a python array.            */
int
downsample(double *input,      /* Input array                               */
           PyArrayObject *out, /* Output array                              */
           int n,              /* Number of elements in the input array     */
           int scale){         /* Resampling factor                         */
  /* - If the scaling factor (f) is an odd number, this routine simply
       performs an averaging of the f points around the resampled value.
     - If f is even, then the boundary points are weighted by one half and
       distributed into the two resampled points.
     - The x coordinates of the first and last values of the input and output
       arrays will coincide.

     The integral area under the curves (the input and output arrays) is
     nearly conserved.

     For example, if the input array is:
       I = {I0 I1 I2 I3 I4 I5 I6}
     The output for a scaling factor f=3 is:
       O0 = (     I0 + I1) / [0.5(f+1)]
       O1 = (I2 + I3 + I4) / [    f   ]
       O2 = (I5 + I6     ) / [0.5(f+1)]
     The output for a scaling factor f=2 is:
       O0 = (         I0 + 0.5 I1) / [0.5(f+1)]
       O1 = (0.5 I1 + I2 + 0.5 I3) / [    f+1 ]
       O2 = (0.5 I3 + I4 + 0.5 I5) / [    f+1 ]
       O3 = (0.5 I5 + I6         ) / [0.5(f+1)]                             */

  int i, j;   /* for-loop indices */
  /* Number of points in the downsampled array:                             */
  int m = 1 + (n-1)/scale;
  //fprintf(stderr, "m=%i\n", m);
  /* Kernel size:                                                           */
  int ks = 2*(scale/2) + 1;

  /* Odd/even flag:                                                         */
  int even = 1;
  if (scale % 2 != 0)
    even = 0;

  /* First point:                                                           */
  INDd(out,0) = 0.0;
  for (i=0; i<ks/2+1; i++)
    INDd(out,0) += input[i];
  if (even == 1)
    INDd(out,0) -= 0.5*input[ks/2];
  INDd(out,0) /= 0.5*(scale+1);

  /* Down-sampling:                                                         */
  for (j=1; j<m-1; j++){
    INDd(out,j) = 0.0;
    for (i=-ks/2; i < ks/2+1; i++){
      INDd(out,j) += input[scale*j + i];
    }
    if (even == 1)
      INDd(out,j) -= 0.5*(input[scale*j-ks/2] + input[scale*j+ks/2]);
    INDd(out,j) /= scale;
  }

  /* Last point:                                                            */
  INDd(out,(m-1)) = 0.0;
  for (i=n-1-ks/2; i<n; i++)
    INDd(out,(m-1)) += input[i];
  if (even == 1)
    INDd(out,(m-1)) -= 0.5*input[n-ks/2];
  INDd(out,(m-1)) /= 0.5*(scale+1);

  return 0;
}

