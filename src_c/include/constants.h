// Copyright (c) 2021 Patricio Cubillos
// Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

/* Math constants:                                                          */
#define PI         (3.141592653589793)       /* pi                          */
#define SQRTLN2    (0.83255461115769775635)  /* sqrt(ln(2))                 */
#define TWOOSQRTPI (1.12837916709551257389)  /* 2*sqrt(pi)                  */
#define SQRTLN2PI  (0.46971863934982566689)  /* sqrt(ln(2*pi))              */

/* Physical constants:                                                      */
#define LS  (2.99792458e10)           /* Light Speed (cm / s)               */
#define KB  (1.380658e-16)            /* Boltzmann constant (erg / K)       */
#define AMU (1.66053886e-24)          /* Atomic Mass unit (g)               */
#define H   (6.6260755e-27)           /* Planck's constant (erg * s)        */
#define EC  (4.8032068e-10)           /* Electronic charge (statcoulomb)    */
#define ME  (9.1093897e-28)           /* Electron mass (g)                  */
#define ATM (1010000.0)       // 1 atm in barye

/* Other constants:                                                         */
#define SIGCTE  (PI*EC*EC/LS/LS/ME)
#define EXPCTE  (H*LS/KB)

/* Constants for Doppler and Lorentz width calculation:                     */
#define FDOP (3.581175136e-7)            /* sqrt(2*KB/AMU) * sqrt(ln2) / c  */
#define FLOR (1.461466451e17)            /* sqrt(2*KB/pi/AMU) / (AMU*c)     */

/* Constants for Voigt-profile calculation:                                 */
#define A1  0.46131350
#define A2  0.19016350
#define A3  0.09999216
#define A4  1.78449270
#define A5  0.002883894
#define A6  5.52534370

#define B1  0.51242424
#define B2  0.27525510
#define B3  0.05176536
#define B4  2.72474500

#define C2  1.4387768775039338
#define C3  8.852821681767784e-13

