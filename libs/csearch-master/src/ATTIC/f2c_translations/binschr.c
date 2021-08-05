/* binschr.f -- translated by f2c (version of 23 April 1993  18:34:30).
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
    shortint ibndp[6250], iangp[9150], itorp[3600], iimpp[3250], ihbp[4100];
} getpar_;

#define getpar_1 getpar_

struct {
    real sc_bond_bld__[200], sc_angle_bld__[200], sc_tors_bld__[200], 
	    sc_offset__[200];
    integer nscatm, nsclmp, nscres, sc_code_bld__[200], sc_atom_part__[76], 
	    sc_clump_part__[31], sc_special__[30], sc_symmetry__[75], 
	    sc_ante1_bld__[200], sc_ante2_bld__[200], sc_ante3_bld__[200], 
	    sc_resname__[30], sc_free_atom__[75], sc_atom_bld__[200];
} sidetop_;

#define sidetop_1 sidetop_

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int binschr_(a, n, key, istart, istop)
real *a;
integer *n;
real *key;
integer *istart, *istop;
{

    static integer ind;


/*     Performs a binary search on the array A for the the KEY. ISTART and
 */
/*    ISTOP will return the array elements which bound KEY. ISTART will be
*/
/*     zero if KEY precedes the first element; ISTOP will be greater than 
*/
/*     N if KEY comes after the last element. */
/*    Moved function parameter definitions to start of function ACRM 07.01
.94*/
/*    Moved before includes as sidetop.inc has a DATA statement ACRM 07.01
.94*/
/* -- Temp */
/* #include "params.inc" */
/* #include "coords.inc" */
/* #include "getpar.inc" */
/* #include "sidetop.inc" */
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
/* *      NAME: GETPAR                                                   *
 */
/* *  FUNCTION: To point to parameters for bond etc energy calculations  *
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

/*  Variable name   Array bounds             Description */
/*  -------------   ------------             ----------- */

/*      IBNDP        MAXBND      Location of bond parameters */
/*      IANGP        MAXANG      Location of angle parameters */
/*      ITORP        MAXTOR      Location of proper torsion parameters */
/*      IIMPP        MAXIMP      Location of improper torsion parameters 
*/
/*      IHBP         MAXHB       Location of hydrogen-bond parameters */

/* -- End Temp */

/* ***********************************************************************
 */
/* *      NAME: SIDETOP                                                  *
 */
/* *  FUNCTION: To declare variables relating to sidechain topology      *
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
/* * 07/01/94   ACRM    Corrected DATA statement to end of file          *
 */
/* ***********************************************************************
 */

/*   Variable     Array bounds            Description */
/*   --------     ------------            ----------- */
/*   IADD              -       Flag for additions to free torsion */
/*   IFIXED            -       Flag indicating fixed torsion */
/*   IFREE             -       Torsion flag indicating independent atom */
/*   NSCATM            -       The number of sidechain atoms */
/*   NSCLMP            -       The number of sidechain clumps */
/*   NSCRES            -       The number of sidechain residues */
/*   SC_ATOM_BLD    MXSCAT     IUPAC name of atom to be constructed */
/*   SC_ANTE1_BLD   MXSCAT     Neighbour IUPAC named atom for building */
/*   SC_ANTE2_BLD   MXSCAT     Next adjacent IUPAC atom for building */
/*   SC_ANTE3_BLD   MXSCAT     Third adjacent IUPAC atom for building */
/*   SC_BOND_BLD    MXSCAT     Bond for building an atom */
/*   SC_ANGLE_BLD   MXSCAT     Construction angle for each atom */
/*  SC_CODE_BLD    MXSCAT     Flag describing the use of the torsion angle
*/
/*  SC_TORS_BLD    MXSCAT     Torsion angle for building (FREE=> dont use)
*/
/*   SC_ATOM_PART   MXCLMP     Partition of atoms into clumps */
/*   SC_CLUMP_PART  MXSRES     Partition of clumps into residues */
/*   SC_FREE_ATOM   MXCLMP     The free atom in the clump */
/*   SC_OFFSET      MXSCAT     Offset angle for torsion ADD procedure */
/*   SC_RESNAME     MXSRES     The residue names */
/*   SC_SPECIAL     MXSRES     TRUE if consistency checks are ignored */
/*   SC_SYMMETRY    MXCLMP     The Rotational symmetry of each clump */
/* Where: */

/*   MXSRES                 Maximum number of sidechain residues */
/*   MXCLMP                 Maximum number of sidechain clumps */
/*   MXSCAT                 Maximum number of sidechain atoms */
    /* Parameter adjustments */
    --a;

    /* Function Body */
    if (*n == 0) {
	*istart = 0;
	*istop = 1;
    } else if (*key < a[1]) {
	*istart = 0;
	*istop = 1;
    } else if (*key > a[*n]) {
	*istart = *n;
	*istop = *n + 1;
    } else {
	*istart = 1;
	*istop = *n;
L50:
	if (*istop - *istart > 1) {
	    ind = (*istop + *istart) / 2;
	    if (*key <= a[ind]) {
		*istop = ind;
	    } else {
		*istart = ind;
	    }
	    goto L50;
	}
    }
    return 0;
} /* binschr_ */

