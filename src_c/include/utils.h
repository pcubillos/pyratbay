// Copyright (c) 2016 Patricio Cubillos and contributors.
// Pyrat Bay is currently proprietary software (see LICENSE).

int
pyramidsearch(PyArrayObject *array, double value, int lo, int hi){
  int scale=1,
      imax;
  /* Out of bounds:                                                         */
  if (value < INDd(array,lo))
    return lo;
  if (INDd(array,hi) < value)
    return hi;

  /* Search up the pyramid:                                                 */
  imax = (lo+scale < hi) ? lo+scale:hi;
  while (INDd(array,imax) < value){
    scale *= 2;
    imax = (lo+scale < hi) ? lo+scale:hi;
  }

  /* Now, binsearch down:                                                   */
  while (imax-lo > 1){
    if (INDd(array,((imax+lo)/2)) < value)
      lo = (imax+lo)/2;
    else
      imax = (imax+lo)/2;
  }
  /* Return indext of closest value:                                        */
  if (fabs(INDd(array,imax)-value) < fabs(INDd(array,lo)-value))
    return imax;
  return lo;
}


int
binsearchapprox(PyArrayObject *array, double value, int lo, int hi){
  /* Last case, value limited between consecutive indices of array:         */
  if (hi-lo == 1){
    /* Return closest array index to value:                                 */
    if (fabs(INDd(array,hi)-value) < fabs(INDd(array,lo)-value))
      return hi;
    return lo;
  }
  /* Compare to middle point and search in corresponding sub-array:         */
  else if (INDd(array,((hi+lo)/2)) > value)
    return binsearchapprox(array, value, lo, (hi+lo)/2);
  else
    return binsearchapprox(array, value, (hi+lo)/2, hi);
}


/* Check if value is in array.  Return the index where it was found, else
   return -1.                                                               */
int
valueinarray(PyArrayObject *array, int value, int arraylen){
  int i;
  for (i=0; i<arraylen; i++){
    if (INDi(array,i) == value)
      return i;
  }
  return -1;
}

/* Find the index of the value with the highest value:                      */
int
imax(PyArrayObject *array){
  int i, max=0;
  for (i=0; i<(int)PyArray_DIM(array,0); i++){
    if (INDi(array, i) > INDi(array, max)){
      max = i;
    }
  }
  return max;
}

/* Downsample an array by an integer factor into a python array.            */
int
downsample(double **input,     /* Input array                               */
           PyArrayObject *out, /* Output array                              */
           int n,              /* Number of elements in the input array     */
           int scale,          /* Resampling factor                         */
           int index){
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
       O1 = (0.5 I1 + I2 + 0.5 I3) / [    f   ]
       O2 = (0.5 I3 + I4 + 0.5 I5) / [    f   ]
       O3 = (0.5 I5 + I6         ) / [0.5(f+1)]                             */

  int i, j;                 /* Auxilliary for-loop indices                  */
  int m = 1 + (n-1)/scale;  /* Number of points in the downsampled array    */
  int ks = 2*(scale/2) + 1; /* Kernel size                                  */
  int even = 1;             /* Odd/even flag:                               */
  if (scale % 2 != 0)
    even = 0;

  /* First point:                                                           */
  for (i=0; i<ks/2+1; i++)
    IND2d(out,index,0) += input[index][i] / (0.5*(scale+1));
  if (even == 1)
    IND2d(out,index,0) -= input[index][ks/2] / (scale+1.0);

  /* Down-sampling:                                                         */
  for (j=1; j<m-1; j++){
    for (i=-ks/2; i < ks/2+1; i++){
      IND2d(out,index,j) += input[index][scale*j + i] / scale;
    }
    if (even == 1)
      IND2d(out,index,j) -= 0.5*(input[index][scale*j-ks/2] +
                                 input[index][scale*j+ks/2])/scale;
  }

  /* Last point:                                                            */
  for (i=n-1-ks/2; i<n; i++)
    IND2d(out,index,(m-1)) += input[index][i] / (0.5*(scale+1));
  if (even == 1)
    IND2d(out,index,(m-1)) -= input[index][n-ks/2] / (scale+1.0);

  return 0;
}


/* Conditional dual (screen/string) message printing.                      */
int
msg(int verb, char *buffer, char *message, ...){
  if (verb < 0)
    return 0;

  va_list args;
  va_start(args, message);

  /* Append formatted text to the end the string:                          */
  vsprintf(buffer+(int)strlen(buffer), message, args);
  /* Print to screen:                                                      */
  vprintf(message, args);
  va_end(args);
  return 0;
}

