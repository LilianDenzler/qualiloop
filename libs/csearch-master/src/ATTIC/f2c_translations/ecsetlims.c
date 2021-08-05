/* ecsetlims.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer dbg_cgen__, dbg_clschn__, dbg_alloc__, dbg_allstk__, dbg_allhp__;
} dbg_;

#define dbg_1 dbg_

struct {
    integer grid_space_grid__, ngridx, ngridy, ngridz;
    real xmn, ymn, zmn, xmx, ymx, zmx, spgridsz, recipgrid;
    integer grid_excluded__, grid_cntnbx__, grid_ingrid__, maxnbx, 
	    grid_nbxa__, freehd, freecls, grid_nexthd__, grid_clshd__, 
	    grid_clstl__, grid_nextcls__, grid_clsatm__, grid_donp__, 
	    grid_accp__, grid_resbya__, grid_qside__, grid_radius__;
    real maxradius;
} grid_;

#define grid_1 grid_

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int ecsetlims_(oc, ecx, avoidxcenter, avoidextent, avoidrmax,
	 avoidvdwrmax, startx, starty, startz, lastx, lasty, lastz, rmax)
real *oc, *ecx, *avoidxcenter, *avoidextent, *avoidrmax, *avoidvdwrmax;
integer *startx, *starty, *startz, *lastx, *lasty, *lastz;
real *rmax;
{
    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer irad, ix, iy, iz;
    static real center[3];
    extern integer chmceil_();


/*     Sets the searching limits for the clump whose x axis is given by */
/*     OC and ECX and whose size is given by the AVOID parameters. RMAX */
/*     is returned as the maximum possible distance. */

/*      IMPLICIT INTEGER(A-Z) */
/* OM OMUPD BNJ 2/9/91 */
/* OM */
/* ACRM added INT definition of CHMCEIL */
/* OM OMUPD BNJ */

/* #include "dbg.inc" */
/* #include "grid.inc" */
/* ***********************************************************************
 */
/* *      NAME: DBG                                                      *
 */
/* *  FUNCTION: To declare the DeBuG level parameters for CONGEN         *
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
/* Declarations */

/*     Data structure for space grid. Pointer variables are prefixed */
/*     by GRID_ in order to avoid naming conflicts. */



    /* Parameter adjustments */
    --ecx;
    --oc;

    /* Function Body */
    center[0] = oc[1] + ecx[1] * *avoidxcenter;
    center[1] = oc[2] + ecx[2] * *avoidxcenter;
    center[2] = oc[3] + ecx[3] * *avoidxcenter;
/* Computing 2nd power */
    r__1 = *avoidrmax;
/* Computing 2nd power */
    r__2 = *avoidextent;
    *rmax = sqrt(r__1 * r__1 + r__2 * r__2) + grid_1.maxradius + *
	    avoidvdwrmax;
    r__1 = *rmax * grid_1.recipgrid;
    irad = chmceil_(&r__1);
    ix = (center[0] - grid_1.xmn) * grid_1.recipgrid + 1;
    iy = (center[1] - grid_1.ymn) * grid_1.recipgrid + 1;
    iz = (center[2] - grid_1.zmn) * grid_1.recipgrid + 1;
    *startx = ix - irad;
    *starty = iy - irad;
    *startz = iz - irad;
    if (*startx < 1) {
	*startx = 1;
    }
    if (*starty < 1) {
	*starty = 1;
    }
    if (*startz < 1) {
	*startz = 1;
    }
    *lastx = ix + irad;
    *lasty = iy + irad;
    *lastz = iz + irad;
    if (*lastx > grid_1.ngridx) {
	*lastx = grid_1.ngridx;
    }
    if (*lasty > grid_1.ngridy) {
	*lasty = grid_1.ngridy;
    }
    if (*lastz > grid_1.ngridz) {
	*lastz = grid_1.ngridz;
    }
    return 0;
} /* ecsetlims_ */

