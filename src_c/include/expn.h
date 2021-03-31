// Copyright (c) 2021 Patricio Cubillos
// Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

#include <math.h>

#define EUL 0.57721566490153286060
#define BIG 1.44115188075855872E+17

#define MACHEP 1.38777878078144567553E-17  /* 2**-56                        */
#define MAXLOG 8.8029691931113054295988E1  /* log(2**127)                   */

/* Exponential integral.
   Taken from:
   https://github.com/scipy/scipy/blob/master/scipy/special/cephes/expn.c   */
static double expn(int n, double x){
  double ans=0.0, r, t, yk, xk;
  double pk, pkm1, pkm2, qk, qkm1, qkm2;
  double psi, z;
  int i, k;
  static double big = BIG;

  if (n < 0){
    printf("Domain error, n (%d) must be > 0.\n", n);
    goto done;
  }
  if (x < 0){
    printf("Domain error, x (%g) must be >= 0.\n", x);
    ans = 1.0/0.0;   /* Return inf                                          */
    goto done;
  }

  if (x > MAXLOG){
    goto done;
  }

  if (x == 0.0) {
    if (n < 2) {
      printf("Domain error, x must be > 0 for n < 2.\n");
      ans = 1.0/0.0;
      goto done;
    }
    else{
      ans = 1.0 / (n - 1.0);
      goto done;
    }
  }

  if (n == 0){
    ans = exp(-x) / x;
    goto done;
  }

  /*  Expansion for large n                                                 */
  if (n > 5000) {
    xk = x + n;
    yk = 1.0 / (xk * xk);
    t = n;
    ans = yk * t * (6.0 * x * x - 8.0 * t * x + t * t);
    ans = yk * (ans + t * (t - 2.0 * x));
    ans = yk * (ans + t);
    ans = (ans + 1.0) * exp(-x) / xk;
    goto done;
  }

  if (x > 1.0)
    goto cfrac;

  /*   Power series expansion                                               */
  psi = -EUL - log(x);
  for (i = 1; i < n; i++)
    psi = psi + 1.0 / i;

  z = -x;
  xk = 0.0;
  yk = 1.0;
  pk = 1.0 - n;
  if (n == 1)
    ans = 0.0;
  else
    ans = 1.0 / pk;
  do {
    xk += 1.0;
    yk *= z / xk;
    pk += 1.0;
    if (pk != 0.0) {
        ans += yk / pk;
    }
    if (ans != 0.0)
      t = fabs(yk / ans);
    else
      t = 1.0;
  }
  while (t > MACHEP);
  k = xk;
  t = n;
  r = n - 1;
  ans = (pow(z, r) * psi / tgamma(t)) - ans;
  goto done;

  /* Continued fraction:                                                    */
  cfrac:
  k = 1;
  pkm2 = 1.0;
  qkm2 = x;
  pkm1 = 1.0;
  qkm1 = x + n;
  ans = pkm1 / qkm1;

  do {
      k += 1;
      if (k & 1) {
          yk = 1.0;
          xk = n + (k - 1) / 2;
      }
      else {
          yk = x;
          xk = k / 2;
      }
      pk = pkm1 * yk + pkm2 * xk;
      qk = qkm1 * yk + qkm2 * xk;
      if (qk != 0) {
          r = pk / qk;
          t = fabs((ans - r) / r);
          ans = r;
      }
      else
          t = 1.0;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      if (fabs(pk) > big) {
          pkm2 /= big;
          pkm1 /= big;
          qkm2 /= big;
          qkm1 /= big;
      }
  }
  while (t > MACHEP);

  ans *= exp(-x);

  done:
    return ans;
}


