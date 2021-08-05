/* abmerfi.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int abmerfi_(p, x, ier)
real *p, *x;
integer *ier;
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    static logical start;
    static real hi;
    extern doublereal abmerf_();
    static real abmerfc, old, bot, cut, top;

/*                **************** */
/*     Calculates the inverse Error Function by doing a binary search */
/*     of calls to the error function. It finds the top and bottom values 
*/
/*     which give the same value of P and takes the average of the two. */
/*     Andrew C.R. Martin     LMB, Oxford. */
/*     This code may be copied and used by anyone providing this notice */
/*     is retained */
/*     PARAMETER (SMALL = 0.3E-7) */
    if (*p == (float)0.) {
	*x = (float)0.;
	*ier = 0;
	return 0;
    }
    if (*p >= (float)1. || *p <= (float)-1.) {
	*ier = 129;
	*x = (float)9.9e24;
	return 0;
    }
    if (*p < (float)0.) {
	*p = -(doublereal)(*p);
    }
/*    First search for the maximum value of X which gives this value for P
*/
    bot = (float)0.;
    top = (float)5.;
    start = TRUE_;
L100:
    cut = bot + (top - bot) / 2;
    abmerfc = abmerf_(&cut);
    if (start) {
	start = FALSE_;
    } else if (cut == old) {
	goto L200;
    }
    old = cut;
    if (abmerfc <= *p) {
	r__1 = cut + (float)3e-21;
	if (abmerf_(&r__1) > *p) {
	    goto L200;
	}
	bot = cut;
    } else {
	top = cut;
    }
    goto L100;
L200:
    hi = cut;
/*     Now find the bottom value satisfying.... */
    top = hi;
    bot = hi - (float)2.;
    start = TRUE_;
    if (bot < (float)0.) {
	bot = (float)0.;
    }
L300:
    cut = bot + (top - bot) / 2;
    if (start) {
	start = FALSE_;
    } else if (cut == old) {
	goto L400;
    }
    old = cut;
    if (cut == hi) {
	*x = hi;
	*ier = 0;
	return 0;
    }
    abmerfc = abmerf_(&cut);
    if (abmerfc >= *p) {
	r__1 = cut - (float)3e-21;
	if (abmerf_(&r__1) < *p) {
	    goto L400;
	}
	top = cut;
    } else {
	bot = cut;
    }
    goto L300;
L400:
    *x = cut + (hi - cut) / 2;
    *ier = 0;
    return 0;
} /* abmerfi_ */

