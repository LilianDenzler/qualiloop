/* gcf.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__4 = 4;
static real c_b30 = (float)3.0000000000000005e-10;
static real c_b59 = (float)3e-20;

/* Subroutine */ int gcf_(gammcf, a, x, gln)
real *gammcf, *a, *x, *gln;
{
    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();
    double log(), exp();

    /* Local variables */
    static real gold, test, g;
    static integer n;
    static real a0, a1, b0, b1, an;
    extern doublereal gammln_();
    static real fac, ana;
    extern /* Subroutine */ int die_();
    static real anf;

    /* Fortran I/O blocks */
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };



/* OM OMUPD BNJ 2/9/91 */




/* OM */

/*     Returns the incomplete gamma function Q(a,x) as a continued */
/*     fraction GAMMCP. */
/*     Also returns ln gamma(a) as GLN */


    *gln = gammln_(a);
    gold = (float)0.;
    a0 = (float)1.;
    a1 = *x;
    b0 = (float)0.;
    b1 = (float)1.;
    fac = (float)1.;
    for (n = 1; n <= 10000; ++n) {
	an = (real) n;
	ana = an - *a;
	a0 = (a1 + a0 * ana) * fac;
	b0 = (b1 + b0 * ana) * fac;
	anf = an * fac;
	a1 = *x * a0 + anf * a1;
	b1 = *x * b0 + anf * b1;
	if (a1 != (float)0.) {
	    fac = (float)1. / a1;
	    g = b1 * fac;
	    if ((float)3e-20 > (r__1 = (g - gold) / g, dabs(r__1))) {
		goto L200;
	    }
	    gold = g;
	}
/* L100: */
    }
    if ((r__1 = (g - gold) / g, dabs(r__1)) < (float)3.0000000000000005e-10) {
	s_wsle(&io___12);
	do_lio(&c__9, &c__1, "Warning in GCF: MAY not have converged", 38L);
	e_wsle();
	s_wsle(&io___13);
	do_lio(&c__9, &c__1, "G =", 3L);
	do_lio(&c__4, &c__1, (char *)&g, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___14);
	do_lio(&c__9, &c__1, "GOLD = ", 7L);
	do_lio(&c__4, &c__1, (char *)&gold, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___15);
	do_lio(&c__9, &c__1, "(G-GOLD)/G =", 12L);
	r__1 = (g - gold) / g;
	do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___16);
	do_lio(&c__9, &c__1, "ABS((G-GOLD)/G) =", 17L);
	r__2 = (r__1 = g - gold, dabs(r__1)) / g;
	do_lio(&c__4, &c__1, (char *)&r__2, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___17);
	do_lio(&c__9, &c__1, "EPS*1.0E10 =", 12L);
	do_lio(&c__4, &c__1, (char *)&c_b30, (ftnlen)sizeof(real));
	e_wsle();
    } else {
	s_wsle(&io___18);
	do_lio(&c__9, &c__1, "Error in GCF: A too large, ITMAX too small", 
		42L);
	e_wsle();
	s_wsle(&io___19);
	do_lio(&c__9, &c__1, "G =", 3L);
	do_lio(&c__4, &c__1, (char *)&g, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___20);
	do_lio(&c__9, &c__1, "GOLD = ", 7L);
	do_lio(&c__4, &c__1, (char *)&gold, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___21);
	do_lio(&c__9, &c__1, "(G-GOLD)/G =", 12L);
	r__1 = (g - gold) / g;
	do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___22);
	do_lio(&c__9, &c__1, "ABS((G-GOLD)/G) =", 17L);
	r__2 = (r__1 = g - gold, dabs(r__1)) / g;
	do_lio(&c__4, &c__1, (char *)&r__2, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___23);
	do_lio(&c__9, &c__1, "EPS =", 5L);
	do_lio(&c__4, &c__1, (char *)&c_b59, (ftnlen)sizeof(real));
	e_wsle();
	die_();
    }

L200:

/* -- Prevent floating point underflows if bounds checking is enabled. */
/*    If TEST is less than '-80.' round down to 0.0 . */

    test = -(doublereal)(*x) + *a * log(*x) - *gln;
    if (test > (float)-80.) {
	*gammcf = exp(test) * g;
    } else {
	*gammcf = (float)0.;
    }

    return 0;
} /* gcf_ */

