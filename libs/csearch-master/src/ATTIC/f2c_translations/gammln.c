/* gammln.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal gammln_(xx)
doublereal *xx;
{
    /* Initialized data */

    static doublereal cof[6] = { 76.18009173,-86.50532033,24.01409822,
	    -1.231739516,.00120858003,-5.36382e-6 };
    static doublereal stp = 2.50662827465;
    static doublereal half = .5;
    static doublereal one = 1.;
    static doublereal fpf = 5.5;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log();

    /* Local variables */
    static integer j;
    static doublereal x, ser, tmp;

/*              ********** */
/*     Returns ln(gamma(xx)) for xx>0. Full accuracy at xx>1. */
/* OM OMUPD BNJ 2/9/91 */
/* OM */
    x = *xx - one;
    tmp = x + fpf;
    tmp = (x + half) * log(tmp) - tmp;
    ser = one;
    for (j = 1; j <= 6; ++j) {
	x += one;
	ser += cof[j - 1] / x;
/* L100: */
    }
    ret_val = tmp + log(stp * ser);
    return ret_val;
} /* gammln_ */

