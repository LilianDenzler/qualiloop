/* ecsetxyzsys.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    real xcart[6150], ycart[6150], zcart[6150], xwork[6150], ywork[6150], 
	    zwork[6150];
} coords_;

#define coords_1 coords_

struct {
    integer dbg_cgen__, dbg_clschn__, dbg_alloc__, dbg_allstk__, dbg_allhp__;
} dbg_;

#define dbg_1 dbg_

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int ecsetxyzsys_(ante, ori, ex, ey, ez)
integer *ante;
real *ori, *ex, *ey, *ez;
{
    /* System generated locals */
    real r__1, r__2, r__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static real dotp, rx, ry;


/*     Constructs a coordinate system for van der Waals repulsion */
/*     avoidance. The origin is on atom ante(2), the x axis is on the */
/*     ante(1)-ante(2) vector, the y axis is the perpendicular component 
*/
/*     of the ante(0)-ante(1) vector, and z axis being the cross product. 
*/

/*      IMPLICIT INTEGER(A-Z) */
/* OM OMUPD BNJ 2/9/91 */
/* OM */

/* #include "params.inc" */
/* #include "coords.inc" */
/* #include "dbg.inc" */
/* ***********************************************************************
 */
/* *      NAME: PARAMS.INC                                               *
 */
/* *  FUNCTION: To declare the main system parameters                    *
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
/*  Parameter Name                    Description */
/*  --------------                    ----------- */
/*     MAXAT      -    Maximum number of atoms */
/*     MAXBND     -    Maximum number of bonds */
/*     MAXANG     -    Maximum number of angles */
/*     MAXTOR     -    Maximum number of proper torsion angles */
/*     MAXIMP     -    Maximum number of improper torsion angles */
/*     MAXHB      -    Maximum number of hydrogen-bonds */
/*     MAXNB      -    Maximum number of non-bond pair exclusions */
/*     MAXRES     -    Maximum number of residues */
/*     MAXSEG     -    Maximum number of segments */
/*     MXDORA     -    Maximum number of hydrogen-bond donors or acceptors
 */
/*     MAXIC      -    Maximum number of internal coordinates */
/*     MAXBP      -    Maximum number of bond parameters */
/*     MAXAP      -    Maximum number of bond angle parameters */
/*     MAXPTP     -    Maximum number of proper torsion parameters */
/*     MAXITP     -    Maximum number of improper torsion parameters */
/*     MAXHBP     -    Maximum number of hydrogen bond parameters */
/*     MAXNBP     -    Maximum number of non-bonded pair parameters */
/*     MAXATU     -    Maximum number of atom types in use */
/*     MAXATT     -    Maximum number of possible atom types */
/*     MXCBUF     -    The main command line buffer size */
/* Declarations */
/* Parameterisations */
/* MAXNBP is related to MAXATU by the formula (MAXATU**2+MAXATU)/2 */
/* ***********************************************************************
 */
/* *      NAME: COORDS.INC                                               *
 */
/* *  FUNCTION: To declare the coordinate arrays                         *
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
/*  Variable name   array bounds           description */
/*  -------------   ------------           ----------- */
/*    XCART           (MAXAT)     X cartesian coordinate for each atom */
/*    YCART           (MAXAT)     Y cartesian coordinate for each atom */
/*    ZCART           (MAXAT)     Z cartesian coordinate for each atom */

/*    XWORK           (MAXAT)     X work space coordinate for each atom */
/*    YWORK           (MAXAT)     Y work space coordinate for each atom */
/*    ZWORK           (MAXAT)     Z work space coordinate for each atom */
/* Declarations */
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

    /* Parameter adjustments */
    --ez;
    --ey;
    --ex;
    --ori;

    /* Function Body */
    ori[1] = coords_1.xcart[ante[2] - 1];
    ori[2] = coords_1.ycart[ante[2] - 1];
    ori[3] = coords_1.zcart[ante[2] - 1];
    ex[1] = coords_1.xcart[ante[2] - 1] - coords_1.xcart[ante[1] - 1];
    ex[2] = coords_1.ycart[ante[2] - 1] - coords_1.ycart[ante[1] - 1];
    ex[3] = coords_1.zcart[ante[2] - 1] - coords_1.zcart[ante[1] - 1];
/* Computing 2nd power */
    r__1 = ex[1];
/* Computing 2nd power */
    r__2 = ex[2];
/* Computing 2nd power */
    r__3 = ex[3];
    rx = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    ex[1] /= rx;
    ex[2] /= rx;
    ex[3] /= rx;
    ey[1] = coords_1.xcart[ante[0] - 1] - coords_1.xcart[ante[1] - 1];
    ey[2] = coords_1.ycart[ante[0] - 1] - coords_1.ycart[ante[1] - 1];
    ey[3] = coords_1.zcart[ante[0] - 1] - coords_1.zcart[ante[1] - 1];
    dotp = ex[1] * ey[1] + ex[2] * ey[2] + ex[3] * ey[3];
    ey[1] -= dotp * ex[1];
    ey[2] -= dotp * ex[2];
    ey[3] -= dotp * ex[3];
/* Computing 2nd power */
    r__1 = ey[1];
/* Computing 2nd power */
    r__2 = ey[2];
/* Computing 2nd power */
    r__3 = ey[3];
    ry = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    ey[1] /= ry;
    ey[2] /= ry;
    ey[3] /= ry;
    ez[1] = ex[2] * ey[3] - ex[3] * ey[2];
    ez[2] = ex[3] * ey[1] - ex[1] * ey[3];
    ez[3] = ex[1] * ey[2] - ex[2] * ey[1];
    return 0;
} /* ecsetxyzsys_ */

