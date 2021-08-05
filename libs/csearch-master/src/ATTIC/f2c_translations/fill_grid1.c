/* fill_grid1.f -- translated by f2c (version of 23 April 1993  18:34:30).
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
    integer natoms, nres, nsegs, nbonds, nangs, nptors, nitors, nhbs, nnbs, 
	    ndonat, naccat, nbpar, napar, nptpar, nitpar, nhbpar, natyps, 
	    nbauto;
} values_;

#define values_1 values_

struct {
    real xcart[6150], ycart[6150], zcart[6150], xwork[6150], ywork[6150], 
	    zwork[6150];
} coords_;

#define coords_1 coords_

struct {
    integer dbg_cgen__, dbg_clschn__, dbg_alloc__, dbg_allstk__, dbg_allhp__;
} dbg_;

#define dbg_1 dbg_

struct {
    real eqbdis[150], bndcon[150], eqang[350], angcon[350], torphs[75], 
	    tormlt[75], torcon[75], eqitan[55], impcon[55], vdwr12[1640], 
	    vdwr6[1640], hbr12[250], hbr10[250], atmpol[100], atneff[100], 
	    vdwrad[100];
    shortint atflag[100];
    integer bndkey[150], angkey[350], impkey[55], torkey[75], nbkey[1640], 
	    hbkey[250], nbcut, dielec, nbflag;
} engpar_;

#define engpar_1 engpar_

struct {
    shortint atbnd1[6250], atbnd2[6250], atang1[9150], atang2[9150], atang3[
	    9150], attor1[3600], attor2[3600], attor3[3600], attor4[3600], 
	    atimp1[3250], atimp2[3250], atimp3[3250], atimp4[3250], hbacpt[
	    1200], hbaan1[1200], hbaan2[1200], hbdonr[1200], hbdhyd[1200], 
	    hbdan1[1200], hbdan2[1200], lstatm[1051], resndx[1050], nbexcl[
	    16150], qmove[6150], atcode[6150];
    integer lstexc[6150], segndx[210]	/* was [10][21] */;
    real segid[20], resid[1050], resnme[1050], atmnme[6150], atchrg[6150], 
	    atmass[6150];
} pstruct_;

#define pstruct_1 pstruct_

/* Table of constant values */

static integer c__0 = 0;


/*     This file contains Fortran subroutines for handling the space grid. */
/*     The additional C code is found in gridc.c. */

/*     Copyright (c) 1987 Robert E. Bruccoleri */
/*     Copying of this software, in whole or in part, is permitted */
/*     provided that the copies are not made for commercial purposes, */
/*     appropriate credit for the use of the software is given, this */
/*     copyright notice appears, and notice is given that the copying */
/*     is by permission of Robert E. Bruccoleri. Any other copying */
/*     requires specific permission. */

/*                     Management of close contacts: */

/*             Close contact searches are performed using Cartesian space */
/*     grid. The grid size is set at a value bigger than any repulsive */
/*     contact whose energy would exceed reasonable values of MAXEVDW. */
/*     The bounds of the grid are determined by the molecular dimensions, */
/*     and the location and size of the loop being searched over. The */
/*     grid consists of pointers to an array of linked lists which store */
/*     the atoms within the cell. The free elements in the array of */
/*     linked list heads and tails are managed by a linked list of free */
/*     elements. A close contact search then requires only searching a */
/*     maximum of 27 cells. It is also relatively quick to remove and add */
/*     atoms to the list which is important as conformations are */
/*     constructed. */

/*             If the space grid overflows, it is expanded to a larger */
/*     size. */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int fill_grid1__(space_grid__, excluded, cntnbx, ingrid, 
	nbxa, nexthd, clshd, clstl, nextcls, clsatm, resbya, radius)
shortint *space_grid__;
logical *excluded;
integer *cntnbx;
logical *ingrid;
shortint *nbxa;
integer *nexthd, *clshd, *clstl, *nextcls, *clsatm, *resbya;
real *radius;
{
    /* System generated locals */
    integer space_grid_dim1, space_grid_dim2, space_grid_offset, nbxa_dim1, 
	    nbxa_offset, i__1, i__2;

    /* Local variables */
    static integer ires;
    extern doublereal maxr_();
    extern /* Subroutine */ int fill2_(), fill4_();
    static integer i, j;
    extern /* Subroutine */ int adatmtgrd_();
    static integer dummy_tail__, ix;
    extern /* Subroutine */ int cprint_(), infrls_(), fixinitgrd_(), cpyprm_()
	    , die_();
    static integer ind, iat;
    extern /* Subroutine */ int filllog_();
    static integer startnb;


/*     Initializes the space grid, along with the free lists used to */
/*     store linked lists of atoms. Also, the excluded array is also */
/*     cleared. */

/* #include "impnone.inc" */
/* #include "grid.inc" */

/*     Data structure for space grid. Pointer variables are prefixed */
/*     by GRID_ in order to avoid naming conflicts. */


/* #include "params.inc" */
/* #include "values.inc" */
/* #include "coords.inc" */
/* #include "dbg.inc" */
/* #include "engpar.inc" */
/* #include "pstruct.inc" */
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
/* ***********************************************************************
 */
/* *      NAME: ENGPAR                                                   *
 */
/* *  FUNCTION: To declare the parameters involved in energy calculations*
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

/* Variable  Array bounds                Description */
/* --------  ------------                ----------- */
/* EQBDIS     MAXBP       List of equilibrium bond lengths */
/* BNDCON     MAXBP       List of bond force constants */
/* EQANG      MAXAP       List of equilibrium bond angles */
/* ANGCON     MAXAP       List of bond angle force constants */
/* TORPHS     MAXPTP      Torsion expression phase shift term */
/* TORMLT     MAXPTP      Torsion angle multiplicity */
/* TORCON     MAXPTP      Torsion angle force constant */
/* EQITAN     MAXITP      Equilibrium improper torsion angle list */
/* IMPCON     MAXITP      Improper torsion angle force constants */
/* VDWR12     MAXNBP      Lennard-Jones R^12 non-bond coefficient */
/* VDWR6      MAXNBP      Lennard-Jones R^6 non-bond coefficient */
/* HBR12      MAXHBP      R to the power 12 H-bond energy terms */
/* HBR10      MAXHBP      R to the power 10 H-bond energy terms */
/* ATMPOL     MAXATT      Atomic polarizability for each atom */
/* ATNEFF     MAXATT      Number of outer shell electrons (effective) */
/* VDWRAD     MAXATT      List of atomic van der Waals radii */
/* ATFLAG     MAXATT      A flag to identify atoms when read from the */
/*                        main topology file */
/* BNDKEY     MAXBP       An array to key a particular bond */
/* ANGKEY     MAXAP       An array to key a particular bond angle */
/* IMPKEY     MAXITP      An array to key a particular improper torsion */
/* TORKEY     MAXPTP      An array to key a particular proper torsion */
/*NBKEY      MAXNBP      An array to key a particular non-bond interaction
*/
/* HBKEY      MAXHBP      An array to key a particular hydrogen-bond */
/* NBCUT        -         The non-bond cutoff distance */
/* DIELEC       -         The dielectric constant */
/* NBFLAG       -         The non-bonded calculation options flag */
/* Declarations */

/* ***********************************************************************
 */
/* *      NAME: PSTRUCT                                                  *
 */
/* *  FUNCTION: To declare the variables used for storing information    *
 */
/* *            related to the whole protein structure. This includes    *
 */
/* *            atomic data and topology data                            *
 */
/* * COPYRIGHT: (C) OXFORD MOLECULAR LIMITED, 1992                       *
 */
/* *---------------------------------------------------------------------*
 */
/* *    AUTHOR: Robert Williams                                          *
 */
/* *      DATE: 06/10/92                                                 *
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
/* Variable  Array Bounds                Description */
/* --------  ------------                ----------- */
/*  ATMNME     MAXAT      List of IUPAC atom names */
/*  ATCODE     MAXAT      List of atom type codes */
/*  ATMASS     MAXAT      List of atomic masses */
/*  ATCHRG     MAXAT      List of atomic charges */
/*  QMOVE      MAXAT      Indicates if atoms can move (0=Y, 1=N) */
/*  ATBND1     MAXBND     List of principal atoms in bonds */
/*  ATBND2     MAXBND     List of second atoms in bonds */
/*  ATANG1     MAXANG     List of first atoms in bond angles */
/*  ATANG2     MAXANG     List of second atoms in bond angles */
/*  ATANG3     MAXANG     List of third atoms in bond angles */
/*  ATTOR1     MAXTOR     List of first atoms in proper torsions */
/*  ATTOR2     MAXTOR     List of second atoms in proper torsions */
/*  ATTOR3     MAXTOR     List of third atoms in proper torsions */
/*  ATTOR4     MAXTOR     List of fourth atoms in proper torsions */
/*  ATIMP1     MAXIMP     List of first atoms in improper torsions */
/*  ATIMP2     MAXIMP     List of second atoms in improper torsions */
/*  ATIMP3     MAXIMP     List of third atoms in improper torsions */
/*  ATIMP4     MAXIMP     List of fourth atoms in improper torsions */
/*  NBEXCL     MAXNB      List of non-bond exclusions, indexed by LSTEXC 
*/
/*  LSTEXC     MAXAT      Gives the last non-bond exclusion in NBEXCL */
/*                        for an atom */
/*  HBACPT     MXDORA     List of hydrogen bond acceptors */
/*  HBAAN1     MXDORA     First antecedent atom to HBACPT */
/*  HBAAN2     MXDORA     Second antecedent atom to HBACPT */
/*  HBDONR     MXDORA     List of heavy atom donors */
/*  HBDHYD     MXDORA     Hydrogen bond donor atom attached hydrogen */
/*  HBDAN1     MXDORA     First antecedant atom to HBDONR */
/*  HBDAN2     MXDORA     Second antecedant atom to HBDONR */
/*  LSTATM   (MAXRES+1)   LSTATM(I+1) gives the last atom in residue I */
/* RESNDX     MAXRES     Index of where residues are in the topology array
s*/
/*  RESNME     MAXRES     List of residue names */
/*  RESID      MAXRES     List of residue indentifications */
/*  SEGID      MAXSEG     List of segment indentifications */
/*  SEGNDX  (10,MAXSEG+1) This array gives information regarding the */
/*                        separation of the structure into segments. */
/*                        The first array argument indicates: */
/*                              1 Total number of residues */
/*                              2 Total number of atoms */
/*                              3 Total number of bonds */
/*                              4 Total number of bond angles */
/*                              5 Total number of proper torsions */
/*                              6 Total number of improper torsions */
/*                              7 Total non-bond exclusions */
/*                              8 Total H-bond donor atoms */
/*                              9 Total H-bond acceptor atoms */
/*                             10 The segment type flag */
/*                        The second argument indicates the last value */
/*                        for the quantity defined by the first argument 
*/
/*                        ie SEGNDX(3,I+1) gives the last bond number */
/*                        for segment I. */

    /* Parameter adjustments */
    --radius;
    --resbya;
    --clsatm;
    --nextcls;
    --clstl;
    --clshd;
    --nexthd;
    nbxa_dim1 = grid_1.maxnbx;
    nbxa_offset = nbxa_dim1 + 1;
    nbxa -= nbxa_offset;
    --ingrid;
    --cntnbx;
    --excluded;
    space_grid_dim1 = grid_1.ngridx;
    space_grid_dim2 = grid_1.ngridy;
    space_grid_offset = space_grid_dim1 * (space_grid_dim2 + 1) + 1;
    space_grid__ -= space_grid_offset;

    /* Function Body */
    i__1 = grid_1.ngridx * grid_1.ngridy * grid_1.ngridz;
    fill2_(&space_grid__[space_grid_offset], &i__1, &c__0);
    i__1 = values_1.natoms;
    for (i = 1; i <= i__1; ++i) {
	excluded[i] = FALSE_;
/* L100: */
    }
    fill4_(&cntnbx[1], &values_1.natoms, &c__0);
    filllog_(&ingrid[1], &values_1.natoms, &c__0);
    i__1 = values_1.natoms;
    for (i = 1; i <= i__1; ++i) {
	if (i == 1) {
	    startnb = 1;
	} else {
	    startnb = pstruct_1.lstexc[i - 2] + 1;
	}
	i__2 = pstruct_1.lstexc[i - 1];
	for (ix = startnb; ix <= i__2; ++ix) {
	    j = pstruct_1.nbexcl[ix - 1];
	    if (j >= 0) {
		if (i < j) {
		    ++cntnbx[i];
		    if (cntnbx[i] > grid_1.maxnbx) {
			cprint_("Error in FILL_GRID1 -- MAXNBX exceeded.", 
				39L);
			die_();
		    }
		    nbxa[cntnbx[i] + i * nbxa_dim1] = j;
		    ++cntnbx[j];
		    if (cntnbx[j] > grid_1.maxnbx) {
			cprint_("Error in FILL_GRID1 -- MAXNBX exceeded.", 
				39L);
			die_();
		    }
		    nbxa[cntnbx[j] + j * nbxa_dim1] = i;
		}
	    }
/* L150: */
	}
/* L200: */
    }
    infrls_(&grid_1.freehd, &dummy_tail__, &nexthd[1], &values_1.natoms);
    infrls_(&grid_1.freecls, &dummy_tail__, &nextcls[1], &values_1.natoms);
    i__1 = values_1.natoms;
    for (ind = 1; ind <= i__1; ++ind) {
	if (coords_1.xcart[ind - 1] != 9999.) {
/* Modified for f2c ACRM 07.01.94 */
/*            CALL ADATMTGRD(%VAL(IND)) */
	    adatmtgrd_(ind);
	}
/* L300: */
    }
    fixinitgrd_();
    i__1 = values_1.nres;
    for (ires = 1; ires <= i__1; ++ires) {
	i__2 = pstruct_1.lstatm[ires];
	for (iat = pstruct_1.lstatm[ires - 1] + 1; iat <= i__2; ++iat) {
	    resbya[iat] = ires;
/* L350: */
	}
/* L400: */
    }
    cpyprm_(engpar_1.vdwrad, engpar_1.atflag, pstruct_1.atcode, &radius[1], &
	    values_1.natoms);
    grid_1.maxradius = maxr_(&radius[1], &values_1.natoms);
    return 0;
} /* fill_grid1__ */

