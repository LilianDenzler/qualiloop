/* delatmfgrd1.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

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

struct {
    real xcart[6150], ycart[6150], zcart[6150], xwork[6150], ywork[6150], 
	    zwork[6150];
} coords_;

#define coords_1 coords_

struct {
    integer dbg_cgen__, dbg_clschn__, dbg_alloc__, dbg_allstk__, dbg_allhp__;
} dbg_;

#define dbg_1 dbg_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int delatmfgrd1_(ind, space_grid__, ingrid, nexthd, clshd, 
	clstl, nextcls, clsatm)
integer *ind;
shortint *space_grid__;
logical *ingrid;
integer *nexthd, *clshd, *clstl, *nextcls, *clsatm;
{
    /* Format strings */
    static char fmt_9001[] = "(\002Error in delatmfgrd1 -- Atom\002,i6,\002 \
deleted twice from CLOSE CONTACTS.\002)";
    static char fmt_9002[] = "(\002Error in delatmfgrd1 --\002,\002 Bad indi\
ces for SPACE_GRID = \002,3i12)";
    static char fmt_9003[] = "(\002Error in delatmfgrd1 -- Atom\002,i12,\002\
 not found in SPACE_GRID.\002)";

    /* System generated locals */
    integer space_grid_dim1, space_grid_dim2, space_grid_offset;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle(), s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    static char buffer[512];
    static integer ix, iy, iz, ispace;
    static logical outofbound;
    extern /* Subroutine */ int cprint_(), die_(), compspcgrd_(), schdls_();

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static icilist io___3 = { 0, buffer, 0, fmt_9001, 512, 1 };
    static cilist io___8 = { 0, 6, 0, 0, 0 };
    static icilist io___9 = { 0, buffer, 0, fmt_9002, 512, 1 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static icilist io___12 = { 0, buffer, 0, fmt_9003, 512, 1 };



/*     Deletes atom IND from the grid. */

/* #include "impnone.inc" */
/* #include "grid.inc" */
/*    THIS FILE ALLOWS FOR IMPLICIT NONE TESTING WITHOUT ERRORS */
/*    DUE TO THE FLECS PREPROCESSOR */


/*     Data structure for space grid. Pointer variables are prefixed */
/*     by GRID_ in order to avoid naming conflicts. */


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
    --clsatm;
    --nextcls;
    --clstl;
    --clshd;
    --nexthd;
    --ingrid;
    space_grid_dim1 = grid_1.ngridx;
    space_grid_dim2 = grid_1.ngridy;
    space_grid_offset = space_grid_dim1 * (space_grid_dim2 + 1) + 1;
    space_grid__ -= space_grid_offset;

    /* Function Body */
    if (! ingrid[*ind]) {
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "Using format 9001", 17L);
	e_wsle();
	s_wsfi(&io___3);
	do_fio(&c__1, (char *)&(*ind), (ftnlen)sizeof(integer));
	e_wsfi();
	cprint_(buffer, 512L);
	die_();
    }
    ingrid[*ind] = FALSE_;
    compspcgrd_(ind, &ix, &iy, &iz, &outofbound);
    if (outofbound) {
	s_wsle(&io___8);
	do_lio(&c__9, &c__1, "Using format 9002", 17L);
	e_wsle();
	s_wsfi(&io___9);
	do_fio(&c__1, (char *)&ix, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iy, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iz, (ftnlen)sizeof(integer));
	e_wsfi();
	cprint_(buffer, 512L);
	die_();
    }
    ispace = space_grid__[ix + (iy + iz * space_grid_dim2) * space_grid_dim1];
    if (ispace == 0) {
	s_wsle(&io___11);
	do_lio(&c__9, &c__1, "Using format 9003", 17L);
	e_wsle();
	s_wsfi(&io___12);
	do_fio(&c__1, (char *)&(*ind), (ftnlen)sizeof(integer));
	e_wsfi();
	cprint_(buffer, 512L);
	die_();
    }
    schdls_(&clshd[ispace], &clstl[ispace], &grid_1.freecls, &nextcls[1], &
	    clsatm[1], ind);
    if (clshd[ispace] == 0) {
	space_grid__[ix + (iy + iz * space_grid_dim2) * space_grid_dim1] = 0;
	nexthd[ispace] = grid_1.freehd;
	grid_1.freehd = ispace;
    }
    return 0;
} /* delatmfgrd1_ */

