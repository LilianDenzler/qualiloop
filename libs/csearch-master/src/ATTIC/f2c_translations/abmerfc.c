/* abmerfc.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static real c_b2 = (float).5;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal abmerfc_(x)
real *x;
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    extern doublereal gammp_(), gammq_();

/*              ******* */
/*     Returns the complementary error function erfc(x) */

/* OM OMUPD BNJ */
/* OM */
    if (*x < (float)0.) {
/* Computing 2nd power */
	r__2 = *x;
	r__1 = r__2 * r__2;
	ret_val = gammp_(&c_b2, &r__1) + (float)1.;
    } else {
/* Computing 2nd power */
	r__2 = *x;
	r__1 = r__2 * r__2;
	ret_val = gammq_(&c_b2, &r__1);
    }
    return ret_val;
} /* abmerfc_ */

