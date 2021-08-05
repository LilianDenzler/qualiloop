/* getuv.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__3 = 3;

/* Subroutine */ int getuv_(x, y, z, ind, u, v)
real *x, *y, *z;
integer *ind;
real *u, *v;
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    extern doublereal vlen_();
    static integer i;
    static real r, dotur;
    extern /* Subroutine */ int ascale_();
    static real dx[3];
    extern doublereal dot_();

/*      IMPLICIT INTEGER(A-Z) */
    /* Parameter adjustments */
    --v;
    --u;
    --ind;
    --z;
    --y;
    --x;

    /* Function Body */
    dx[0] = x[ind[2]] - x[ind[1]];
    dx[1] = y[ind[2]] - y[ind[1]];
    dx[2] = z[ind[2]] - z[ind[1]];
    r = vlen_(dx, &c__3);
    r__1 = (float)1. / r;
    ascale_(dx, &r__1, &u[1], &c__3);
    dx[0] = x[ind[3]] - x[ind[2]];
    dx[1] = y[ind[3]] - y[ind[2]];
    dx[2] = z[ind[3]] - z[ind[2]];
    dotur = dot_(&u[1], dx, &c__3);
    for (i = 1; i <= 3; ++i) {
	v[i] = dx[i - 1] - dotur * u[i];
/* L100: */
    }
    r = vlen_(&v[1], &c__3);
    r__1 = (float)1. / r;
    ascale_(&v[1], &r__1, &v[1], &c__3);
    return 0;
} /* getuv_ */

