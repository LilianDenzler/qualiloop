/* cartx2.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer natoms, nres, nsegs, nbonds, nangs, nptors, nitors, nhbs, nnbs, 
	    ndonat, naccat, nbpar, napar, nptpar, nitpar, nhbpar, natyps, 
	    nbauto;
} values_;

#define values_1 values_

/* Table of constant values */

static integer c__1 = 1;


/*     Alternate implementation of CARTX which performs the basic BILDER */
/*     function, the construction of an atom's position from the three */
/*     atoms previous on the chain using the bond length, bond angle, and */
/*     torsion angle. We construct the position of atom I4 from I1, I2, */
/*     and I3. The bond is between I3 and I4, the angle is between I2, */
/*     I3, and I4, and the torsion is through all four. X, Y, Z are the */
/*     coordinate arrays. This implementation is faster than CARTX. */

/*     The basic algorithm is to translate and rotate the atoms so that */
/*     X3 is on the origin, X2 on is the negative x axis, and X1 is on */
/*     the XY plane. The position of X4 is easily generated, and then the */
/*     inverse transformation is performed to bring X4 to its proper */
/*     place. */

/* OM OMUPD BNJ 2/9/91 */
/* Subroutine */ int cartx2_(x, y, z, i1, i2, i3, i4, bond, theta, phi)
real *x, *y, *z;
integer *i1, *i2, *i3, *i4;
real *bond, *theta, *phi;
{
    /* Format strings */
    static char fmt_9001[] = "(\002 Atom 1 index = \002,i5,3f10.3)";
    static char fmt_9002[] = "(\002 Atom 2 index = \002,i5,3f10.3)";
    static char fmt_9003[] = "(\002 Atom 3 index = \002,i5,3f10.3)";

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();
    double sin(), cos(), sqrt();

    /* Local variables */
    static real lxz2, xz2o, cphi, bsin, ctht, sphi, lxz22, stht, l2, x1, x2, 
	    x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, ovlyz1, ovlxz2;
    static char buffer[100];
    extern /* Subroutine */ int cprint_();
    static real y1o, z1o, x2o, z2o, y2o, xx1, xx4, yy4, zz4;
    extern /* Subroutine */ int die_();
    static real ovl2, lyz1;

    /* Fortran I/O blocks */
    static icilist io___2 = { 0, buffer, 0, fmt_9001, 100, 1 };
    static icilist io___3 = { 0, buffer, 0, fmt_9002, 100, 1 };
    static icilist io___4 = { 0, buffer, 0, fmt_9003, 100, 1 };




/* -- Includes:- */

/* #include "values.inc" */

/* -- Declarations:- */

/* ***********************************************************************
 */
/* *      NAME: VALUES                                                   *
 */
/* *  FUNCTION: To declare the most frequently used system values        *
 */
/* * COPYRIGHT: (C) OXFORD MOLECULAR LIMITED, 1992                       *
 */
/* *---------------------------------------------------------------------*
 */
/* *    AUTHOR: Robert Williams                                          *
 */
/* *      DATE: 05/10/92                                                 *
 */
/* *---------------------------------------------------------------------*
 */
/* *    INPUTS:                                                          *
 */
/* *   OUTPUTS:                                                          *
 */
/* *    LOCALS:                                                          *
 */
/* *   GLOBALS:                                                          *
 */
/* *     CALLS:                                                          *
 */
/* *---------------------------------------------------------------------*
 */
/* * MODIFICATION RECORD                                                 *
 */
/* * DD/MM/YY   INITS   COMMENTS                                         *
 */
/* ***********************************************************************
 */
/*    Variable              Description */
/*    --------              ----------- */

/*     NATOMS        The number of atoms in the system */
/*     NRES          The number of residues */
/*     NSEGS         The number of segments */
/*     NBONDS        The number of bonds in the system */
/*     NANGS         The number of bond angles in the system */
/*     NPTORS        The number of proper torsion angles */
/*     NITORS        The number of improper torsion angles */
/*     NHBS          The number of hydrogen bonds */
/*     NNBS          The number of non-bond pairs */
/*     NDONAT        The number of hydrogen bond donor atoms */
/*     NACCAT        The number of hydrogen bond acceptor atoms */
/*     NBPAR         The number of bond parameters */
/*     NAPAR         The number of angle parameters */
/*     NPTPAR        The number of proper torsion parameters */
/*     NITPAR        The number of improper torsion parameters */
/*     NHBPAR        The number of hydrogen parameters */
/*     NATYPS        The number of atom types */
/*     NBAUTO        Flag to generate non-bonded exclusions */

/*     PI            The value of PI! */
/*     RAD120        120 degrees expressed as radians */
/*     DTORAD        Degree TO Radian conversion parameter */
/*     LARGNUM       Largest acceptable real number */
/*     LARGINT       Largest acceptable integer */
/*     ANUM          The coordinates of dummy atoms */
/* Declarations */



/* -- Code:- */

    /* Parameter adjustments */
    --z;
    --y;
    --x;

    /* Function Body */
    if (x[*i1] == 9999. || x[*i2] == 9999. || x[*i3] == 9999.) {
	cprint_("Error in CARTX2 -- One of the antecedent atoms is undefined."
		, 60L);
	s_wsfi(&io___2);
	do_fio(&c__1, (char *)&(*i1), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x[*i1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y[*i1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z[*i1], (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 100L);
	s_wsfi(&io___3);
	do_fio(&c__1, (char *)&(*i2), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x[*i2], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y[*i2], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z[*i2], (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 100L);
	s_wsfi(&io___4);
	do_fio(&c__1, (char *)&(*i3), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x[*i3], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y[*i3], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z[*i3], (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 100L);
	die_();
    }

/*     Generate position of X4 with everything else easily lined up. */

    stht = sin(3.141592653589794 - *theta);
    ctht = cos(3.141592653589794 - *theta);
    sphi = sin(*phi);
    cphi = cos(*phi);
    bsin = *bond * stht;
    x4 = *bond * ctht;
    y4 = bsin * cphi;
    z4 = bsin * sphi;

/*     Translate X1 and X2 so that X3 would be at the origin. */

    x3 = x[*i3];
    y3 = y[*i3];
    z3 = z[*i3];
    x1 = x[*i1] - x3;
    y1 = y[*i1] - y3;
    z1 = z[*i1] - z3;
    x2 = x[*i2] - x3;
    y2 = y[*i2] - y3;
    z2 = z[*i2] - z3;

/*     Rotate X1 by rotation of X2 to the origin. */

    lxz22 = x2 * x2 + z2 * z2;
    l2 = sqrt(lxz22 + y2 * y2);
    lxz2 = sqrt(lxz22);
    if (l2 < (float)1e-7) {
	cprint_("0***** Warning from CARTX ***** Atom 2 and Atom 3 are too c\
lose.", 64L);
	ovl2 = (float)1e7;
    } else {
	ovl2 = (float)1. / l2;
    }
    if (lxz2 < (float)1e-7) {
	xx1 = x1;
	x2o = (float)1.;
	z2o = (float)0.;
    } else {
	ovlxz2 = (float)1. / lxz2;
	x2o = x2 * ovlxz2;
	z2o = z2 * ovlxz2;
	xx1 = x1 * x2o + z1 * z2o;
	z1 = z1 * x2o - x1 * z2o;
    }
    xz2o = lxz2 * ovl2;
    y2o = y2 * ovl2;
    x1 = -(doublereal)xx1 * xz2o - y1 * y2o;
    y1 = xx1 * y2o - y1 * xz2o;

/*     Rotate X4 by inverse of rotation which takes the transformed X1 to 
*/
/*     the XY plane by rotation about the x axis. */

    lyz1 = sqrt(y1 * y1 + z1 * z1);
    ovlyz1 = (float)1. / lyz1;
    y1o = y1 * ovlyz1;
    z1o = z1 * ovlyz1;
    yy4 = y1o * y4 - z1o * z4;
    zz4 = y1o * z4 + z1o * y4;

/*     Rotate X4 by inverse of X2 rotation to -X axis. */

    xx4 = y2o * yy4 - xz2o * x4;
    y4 = -(doublereal)xz2o * yy4 - y2o * x4;
    x4 = x2o * xx4 - z2o * zz4;
    z4 = z2o * xx4 + x2o * zz4;
    x[*i4] = x4 + x3;
    y[*i4] = y4 + y3;
    z[*i4] = z4 + z3;

    return 0;


} /* cartx2_ */

