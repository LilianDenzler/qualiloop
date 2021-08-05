/* gser.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int gser_(gamser, a, x, gln)
real *gamser, *a, *x, *gln;
{
    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();
    double log(), exp();

    /* Local variables */
    static integer n;
    static real ap;
    extern doublereal gammln_();
    extern /* Subroutine */ int die_();
    static real del, sum;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };


/*                ******************** */
/*     Returns the incomplete gamma function P(a,x) as series GAMSER. */
/*     Also returns ln gamma(a) as GLN */
/* OM is this a function call or what? */
    *gln = gammln_(a);
    if (*x <= (float)0.) {
	if (*x < (float)0.) {
	    s_wsle(&io___1);
	    do_lio(&c__9, &c__1, "Error in GSER", 13L);
	    e_wsle();
	    die_();
	}
	*gamser = (float)0.;
	return 0;
    }
    ap = *a;
    sum = (float)1. / *a;
    del = sum;
    for (n = 1; n <= 10000; ++n) {
	ap += (float)1.;
	del = del * *x / ap;
	sum += del;
	if (dabs(del) < dabs(sum) * (float)3e-20) {
	    goto L200;
	}
/* L100: */
    }
    s_wsle(&io___6);
    do_lio(&c__9, &c__1, "Error in GSER: A too large, ITMAX too small", 43L);
    e_wsle();
    die_();
L200:
    *gamser = sum * exp(-(doublereal)(*x) + *a * log(*x) - *gln);
    return 0;
} /* gser_ */

