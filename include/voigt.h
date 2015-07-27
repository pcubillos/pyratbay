/*
 * voigt.c - Functions to compute the Voigt profile
 *
 * Copyright (C) 2004 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

/*
  voigt.c: Functions to return a Voigt profile according to the numerical
  approximation described by Pierlusi et al. in
  J. Quant. Spectrosc. Radiat. Transfer., Vol 18 pp.555

   Normalized line shape is given by:
    \Psi(x,y) = \frac{y}{\pi}
                \int_{-\infty}^{\infty} \frac{\exp(-t^2)}{y^2+(x-t)^2} {\rm d}t
   However, all the functions to be used will be approximations in
   different regions to the function:
     Psi(x,y) = Re[w(z=x+iy)]
              = Re[exp(-z^2)(1-erf(-iz))]

   (c) Patricio Rojo 2003                                                  */

  /* TD: Set iteration limit from fitted function and not convergence
     check!  */
  /* TD: check that function using voigt array hanndles correctly if the
     item happens to be exactly in between two bins: it should be 0.5
     the contribution from its of its boundary bin in the first and
     latter \emph{center shift position}*/
//\omitfh

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "constants.h"

#ifndef _PROFILE_H
#define _PROFILE_H

#define VOIGT_QUICK 0x00001   //Quick integration.
extern int _voigt_maxelements;

#endif /* _PROFILE_H */
/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  */

#define MAXCONV 61
/* ferf[n]= 1/(n!(2n+1)) */
static double ferf[MAXCONV] = {
  1.000000000000000000000,    //0
  0.333333333333333333333,    //1
  0.100000000000000000000,    //2
  2.38095238095238095238e-2,  //3
  4.62962962962962962963e-3,  //4
  7.57575757575757575758e-4,  //5
  1.06837606837606837607e-4,  //6
  1.32275132275132275132e-5,  //7
  1.45891690009337068161e-6,  //8
  1.45038522231504687645e-7,  //9
  1.31225329638028050726e-8,  //10
  1.08922210371485733805e-9,  //11
  8.35070279514723959168e-11, //12
  5.94779401363763503681e-12, //13
  3.95542951645852576340e-13, //14
  2.46682701026445692771e-14, //15
  1.44832646435981372650e-15, //16
  8.03273501241577360914e-17, //17
  4.22140728880708823303e-18, //18
  2.10785519144213582486e-19, //19
  1.00251649349077191670e-20, //20
  4.55184675892820028624e-22, //21
  1.97706475387790517483e-23, //22
  8.23014929921422135684e-25, //23
  3.28926034917575173275e-26, //24
  1.26410789889891635220e-27, //25
  4.67848351551848577373e-29, //26
  1.66976179341737202699e-30, //27
  5.75419164398217177220e-32, //28
  1.91694286210978253077e-33, //29
  6.18030758822279613746e-35, //30
  1.93035720881510785656e-36, //31
  5.84675500746883629630e-38, //32
  1.71885606280178362397e-39, //33
  4.90892396452342296700e-41, //34
  1.36304126177913957635e-42, //35
  3.68249351546114573519e-44, //36
  9.68728023887076175384e-46, //37
  2.48306909745491159104e-47, //38
  6.20565791963739670594e-49, //39
  1.51310794954121709805e-50, //40
  3.60157930981012591661e-52, //41
  8.37341968387228154283e-54, //42
  1.90254122728987952724e-55, //43
  4.22678975419355257584e-57, //44
  9.18642950239868569596e-59, //45
  1.95410258232417110410e-60, //46
  4.07013527785325672298e-62, //47
  8.30461450592911058168e-64, //48
  1.66058051345108993284e-65, //49
  3.25539546201302778914e-67, //50
  6.25918411694871134025e-69, //51
  1.18076183891157008800e-70, //52
  2.18621042295388572103e-72, //53
  3.97425272266506578576e-74, //54
  7.09571739181805357327e-76, //55
  1.24466597738907071213e-77, //56
  2.14564844309633852739e-79, //57
  3.63615636540051474579e-81, //58
  6.05939744697137480783e-83, //59
  9.93207019544894768776e-85
};
int _voigt_maxelements=99999;
int _voigt_computeeach=10;

#define NFCN(x,y) (x<1?15:(int)(6.842*x+8.0))

static inline int
meanintegSimp(double *in,
              double *out,
              int ni,
              int no,
              int ipo,
              double d);


static inline int
meanintegTrap(double *in,
              double *out,
              int ni,
              int no,
              int ipo,
              double d);


static inline void
voigtxy(double x,
        double y,
        double *res,
        double eps,
        double alphaD){
  int i, n;
  long double or, nr, oi, ni, ar, ai; 
  const long double x2y2 = x*x-y*y;
  const long double xy2 = 2*x*y;
  const long double cosxy = cosl(xy2);
  const long double sinxy = sinl(xy2);

#ifdef _VOIGT_USE_CONVERG
  if(eps > 0)
    n = MAXCONV;
  else 
#endif /* Precision convergence */
    n = NFCN(x, y) + 1;

  if(x<3 && y<1.8){ /* Region I */
    ar = or =  y;
    ai = oi = -x;
    i = 1;
    while(1){
      ni = or*xy2  + oi*x2y2;
      nr = or*x2y2 - oi*xy2 ;
#ifdef _VOIGT_USE_CONVERG
      if (n == MAXCONV){
        if (i > n){
          fprintf(stderr, "%s:: No convergence after %i iterations in "
                          "VOIGTXY() calculations (voigt.c file)\n",
                          __FILE__, MAXCONV);
          break;
        }
        if (fabs(ferf[i]*(cosxy*nr + sinxy*ni))<eps)
          break;
      }
      else 
#endif /* Precision convergence */
        if (i > n)
          break;
      ai += ni*ferf[i];
      ar += nr*ferf[i];

      oi = ni;
      or = nr;
      i++;
    }
    (*(res)) = SQRTLN2PI / alphaD * exp(-x2y2) *
               (cosxy * (1-ar*TWOOSQRTPI) - sinxy * ai * TWOOSQRTPI);
  }
  else if(x<5 && y<5){    /* Region II  */
    ar = xy2 * xy2;
    nr = xy2 * x;
    ni = x2y2 - A2;
    ai = x2y2 - A4;
    oi = x2y2 - A6;
    (*(res)) = SQRTLN2PI/alphaD*(A1*((nr-ni*y)/(ni*ni+ar)) +
                                 A3*((nr-ai*y)/(ai*ai+ar)) +
                                 A5*((nr-oi*y)/(oi*oi+ar)) );
  }
  else{                   /* Region III */
    ar = xy2 * xy2;
    nr = xy2 * x;
    ni = x2y2 - B2;
    ai = x2y2 - B4;
    (*(res)) = SQRTLN2PI/alphaD*(B1*((nr-ni*y)/(ni*ni+ar)) +
                                 B3*((nr-ai*y)/(ai*ai+ar)) );
  }
} __attribute__((always_inline))


/* Compute Voigt profile on equispaced grid.
   Return: 1 on success                                                     */
inline int 
voigtn(int nwn,        /* Number of wavenumber bins of the profile          */
       double dwn,     /* Profile half-width (in cm-1)                      */
       double alphaL,  /* Lorentz width                                     */
       double alphaD,  /* Doppler width                                     */
       double *vpro,   /* Array (nwn) where to store the profile            */
       double eps,     /* Precision to which the profile is to be calculated.
                          If negative, do a fixed number of iterations      */
       int flags){     /* VOIGT_QUICK flag to perform a quick integration   */

  /* Initialization:                                                        */
  double *aint; /* Array to put the profile                                 */
  int i, osamp; /* Oversampling factor                                      */
  double y = SQRTLN2 * alphaL / alphaD,     /* Normalized width             */
         x, /* = sqrt(ln 2) * |nu-nu0|/alphaD, Normalized position          */
         ddwn = 2.0 * dwn / (nwn-1);        /* Wavenumber sampling interval */

  /* If the Doppler width is resolved by less than 50 points, calculate the
     profile over an onversampled array, and then bin down.                 */

  int nint = 50; /* Minimum number of points calculated per Doppler width   */
  /* Wavenumber sampling interval assuming nint:                            */
  double dint = alphaD / (nint - 1);

  /* If requested sampling interval is finer than dint, or if VOIGT_QUICK
     is requested, do not oversample:                                       */
  if( ddwn < dint || (flags & VOIGT_QUICK) ){
    osamp = 1;
    dint = ddwn;     /* Calculated profile sampling interval                */
    nint = nwn + 1;  /* Number of profile samples to calculate              */
  }
  /* Else, oversample:                                                      */
  else{
    osamp = (int)(ddwn / dint) + 1;
    if(osamp & 1)
      osamp++;
    /* Enforce nint to an odd number (for Simpson integration):             */
    nint = nwn * osamp + 1;        /* Number of profile samples calculated  */
    dint = 2.0 * dwn / (nint - 1); /* Calculated-profile sampling interval  */
  }

  /* Initialize aint array:                                                 */
  if((aint=(double *)calloc(nint, sizeof(double)))==NULL){
    fprintf(stderr, "\nUnable to allocate memory for Voigt profile.\n");
    exit(EXIT_FAILURE);
  }

  /* Evaluate the Voigt profile point by point:                             */
  for(i=0; i < nint; i++){
    x = SQRTLN2 * fabs(dint*i-dwn) / alphaD;
    voigtxy(x, y, aint+i, eps, alphaD);
  }

  /* Quick integration, return the value at the beginning of each bin:      */
  if(flags & VOIGT_QUICK){ 
    for(i=0; i<nwn; i++)
      vpro[i] = aint[i];
  }

  /* Slow integration, calculate the number of fine-points per regular 
     point, note that is one more than the ratio:                           */
  else{
    /* Is (osamp + 1)  an odd number?
       If so, use Simpson's integration, otherwise use trapezoidal:         */
    if ((osamp+1) & 1)
      meanintegSimp(aint, vpro, nint, nwn, osamp+1, dint);
    else
      meanintegTrap(aint, vpro, nint, nwn, osamp+1, dint);
  }
  free(aint);

  /* On success return 1:                                                   */
  return 1;
}


/* Simpson integration
   Return: 1 on success                                                     */
static inline int
meanintegSimp(double *in,  /* Input array                                   */
              double *out, /* Output array                                  */
              int ni,      /* Number of input elements                      */
              int no,      /* Number of output elements                     */
              int ipo,     /* Number of inputs per output                   */
              double d){   /* bin width                                     */

  /* Indices for the in array. To speed up, \vr{ipo} will be the
     last index in the sub-arrays instead of the number of elements:        */
  int i;
  ipo--;

  for(; no; no--){
    /* Integrate:                                                           */
    *out = 0;
    for(i=1; i<ipo; i+=2)
      *out += in[i];
    *out *= 2;
    for(i=2; i<ipo; i+=2)
      *out += in[i];
    *out *= 2;
    *out += *in + in[ipo];
    *out /= (ipo*3.0);
    /* Advance the input array to the next bin and advance in one the
       return array */
    in += ipo;
    out++;
  }

  return 1;
}


/* Trapezoid integration
   Returns 1 on success.                                                    */
static inline int
meanintegTrap(double *in,  /* Input Array                                   */
              double *out, /* Output Array                                  */
              int ni,      /* Number of input elements                      */
              int no,      /* Number of output elements                     */
              int ipo,     /* Number of inputs per output                   */
              double d){   /* bin width                                     */

  int i;  /* for-loop index                                                 */
  ipo--;  /* Last index in the sub-arrays */

  /* Integrate each bin:                                                    */
  for(; no; no--){
    *out = 0;
    for(i=1; i<ipo; i++)
      *out += in[i];
    *out = (*out + (in[0] + in[ipo])/2.0) / (double)ipo;
    /* Move input/output pointers to the next bin:                          */
    in += ipo;
    out++;
  }

  return 1;
}

#undef SQRTLN2
#undef TWOOSQRTPI
#undef MAXCONV
#undef NFCN
#undef VOIGTXY

#undef A1
#undef A2
#undef A3
#undef A4
#undef A5
#undef A6

#undef B1
#undef B2
#undef B3
#undef B4
