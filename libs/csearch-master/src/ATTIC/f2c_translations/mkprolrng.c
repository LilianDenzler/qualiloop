/* mkprolrng.f -- translated by f2c (version of 23 April 1993  18:34:30).
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

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int mkprolrng_(x, y, z, ind_c__, ind_n__, ind_ca__, ind_cb__,
	 ind_cg__, ind_cd__, phi, procons, proconsphi, nprocons)
real *x, *y, *z;
integer *ind_c__, *ind_n__, *ind_ca__, *ind_cb__, *ind_cg__, *ind_cd__;
real *phi, *procons, *proconsphi;
integer *nprocons;
{
    /* Format strings */
    static char fmt_9001[] = "(\0020Error in MKPROLRNG -- Phi value =\002,1p\
g14.7,\002 is out of bounds.\002)";

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    static real bond[3];
    static integer stop;
    static real theta[3];
    static integer start;
    extern /* Subroutine */ int cartx2_();
    static char buffer[100];
    extern /* Subroutine */ int cprint_(), die_(), binschr_();
    extern doublereal interpa_();
    static real phi1, phi2, torsion[3];

    /* Fortran I/O blocks */
    static icilist io___4 = { 0, buffer, 0, fmt_9001, 100, 1 };



/*     Given a value for the PHI torsion, and pointers to various atoms */
/*     of the proline residue (in IND_C, IND_N, etc.), and the proline */
/*     constructor arrays, the positions of CB, CG and CD are */
/*     constructed. IND_C points to the carbonyl carbon of the previous */
/*     residue. PROCONS contains the internal coordinates for the three */
/*     ring atoms in the order, CA-CB, N-CA-CB, C-N-CA-CB, CB-CG, */
/*     CA-CB-CG, N-CA-CB-CG, CG-CD, CB-CG-CD, CA-CB-CG-CD. PROCONSPHI */
/*     contains the values for PHI to which the constructor information */
/*     applies, and NPROCONS contains the number of Phi's in the tables. 
*/

/* #include "values.inc" */

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
    /* Parameter adjustments */
    --proconsphi;
    procons -= 10;
    --z;
    --y;
    --x;

    /* Function Body */
    binschr_(&proconsphi[1], nprocons, phi, &start, &stop);
    if (start < 1 || stop > *nprocons) {
	s_wsfi(&io___4);
	do_fio(&c__1, (char *)&(*phi), (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 100L);
	die_();
    }
    phi1 = proconsphi[start];
    phi2 = proconsphi[stop];
    bond[0] = interpa_(&phi1, &procons[start * 9 + 1], &phi2, &procons[stop * 
	    9 + 1], phi);
    theta[0] = interpa_(&phi1, &procons[start * 9 + 2], &phi2, &procons[stop *
	     9 + 2], phi);
    torsion[0] = interpa_(&phi1, &procons[start * 9 + 3], &phi2, &procons[
	    stop * 9 + 3], phi);
    bond[1] = interpa_(&phi1, &procons[start * 9 + 4], &phi2, &procons[stop * 
	    9 + 4], phi);
    theta[1] = interpa_(&phi1, &procons[start * 9 + 5], &phi2, &procons[stop *
	     9 + 5], phi);
    torsion[1] = interpa_(&phi1, &procons[start * 9 + 6], &phi2, &procons[
	    stop * 9 + 6], phi);
    bond[2] = interpa_(&phi1, &procons[start * 9 + 7], &phi2, &procons[stop * 
	    9 + 7], phi);
    theta[2] = interpa_(&phi1, &procons[start * 9 + 8], &phi2, &procons[stop *
	     9 + 8], phi);
    torsion[2] = interpa_(&phi1, &procons[start * 9 + 9], &phi2, &procons[
	    stop * 9 + 9], phi);
    cartx2_(&x[1], &y[1], &z[1], ind_c__, ind_n__, ind_ca__, ind_cb__, bond, 
	    theta, torsion);
    cartx2_(&x[1], &y[1], &z[1], ind_n__, ind_ca__, ind_cb__, ind_cg__, &bond[
	    1], &theta[1], &torsion[1]);
    cartx2_(&x[1], &y[1], &z[1], ind_ca__, ind_cb__, ind_cg__, ind_cd__, &
	    bond[2], &theta[2], &torsion[2]);
    return 0;
} /* mkprolrng_ */

