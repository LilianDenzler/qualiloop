/* merfi.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__4 = 4;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int merfi_(p, x, ier)
real *p, *x;
integer *ier;
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static logical start;
    static real hi;
    extern doublereal abmerf_();
    static real abmerfc, old, bot, cut, top;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___2 = { 0, 6, 0, 0, 0 };


/*                ************** */
/*     Calculates the inverse Error Function by doing a binary search */
/*     of calls to the error function. It finds the top and bottom values 
*/
/*     which give the same value of P and takes the average of the two. */
/*     Andrew C.R. Martin     LMB, Oxford. */
/*     This code may be copied and used by anyone providing this notice */
/*     is retained */
    if (*p == (float)0.) {
	*x = (float)0.;
	*ier = 0;
	return 0;
    }
    if (*p > (float)1. || *p < (float)-1.) {
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "Error in MERFI: P out of range", 30L);
	e_wsle();
	s_wsle(&io___2);
	do_lio(&c__9, &c__1, "P = ", 4L);
	do_lio(&c__4, &c__1, (char *)&(*p), (ftnlen)sizeof(real));
	e_wsle();
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
	r__1 = cut + (float)3e-8;
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
	r__1 = cut - (float)3e-8;
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
} /* merfi_ */

