/* gammp.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal gammp_(a, x)
real *a, *x;
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    extern /* Subroutine */ int gser_();
    static real gammcf, gamser;
    extern /* Subroutine */ int gcf_(), die_();
    static real gln;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };


/* OM OMUPD BNJ 2/9/91 */
/* OM */
/*     Returns the incomplete gamma function P(a,x) */
    if (*x < (float)0. || *a <= (float)0.) {
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "Error in GAMMP", 14L);
	e_wsle();
	die_();
    }
    if (*x < *a + (float)1.) {
	gser_(&gamser, a, x, &gln);
	ret_val = gamser;
    } else {
	gcf_(&gammcf, a, x, &gln);
	ret_val = (float)1. - gammcf;
    }
    return ret_val;
} /* gammp_ */

