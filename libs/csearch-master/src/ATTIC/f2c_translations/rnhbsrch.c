/* rnhbsrch.f -- translated by f2c (version of 23 April 1993  18:34:30).
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
    real grid2, cutnb2, cuthb2, epsilon;
    logical cons_die__;
    real cos_cuthba__;
    integer ioff[100], parm_nop__, aamapp, glymapp, promapp, naamap, nglymap, 
	    npromap, proconsp, proconsphip, eproconsp, nprocons;
    char ctitle[800]	/* was [80][10] */;
    integer nctitl, glymapu, alamapu, promapu, proconsu, stunit;
    real glyemax, alaemax, proemax, eringpro;
    logical ignoreevdw;
    integer maxleaf, restart_stlen__, restart_st__;
    logical save_coor__;
    integer nconsp, consp, savex, savey, savez;
    logical sidehits_opt__;
    real maxdt_def__;
} cg_;

#define cg_1 cg_

struct {
    integer dof_headp__;
} dof_head__;

#define dof_head__1 dof_head__

struct {
    integer dof_tailp__;
} dof_tail__;

#define dof_tail__1 dof_tail__

struct {
    integer confnum;
} confnum_;

#define confnum_1 confnum_

struct {
    integer leafnum;
} leafnum_;

#define leafnum_1 leafnum_

struct {
    integer top_level_envp__;
} top_level_env__;

#define top_level_env__1 top_level_env__

struct {
    integer dummy_sidehitsp__;
} dummy_sidehits__;

#define dummy_sidehits__1 dummy_sidehits__

struct {
    real xcart[6150], ycart[6150], zcart[6150], xwork[6150], ywork[6150], 
	    zwork[6150];
} coords_;

#define coords_1 coords_

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

struct {
    integer dbg_cgen__, dbg_clschn__, dbg_alloc__, dbg_allstk__, dbg_allhp__;
} dbg_;

#define dbg_1 dbg_

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int rnhbndsrch_(startx, starty, startz, lastx, lasty, lastz, 
	space_grid__, clshd, clsatm, excluded, search_mode__, donp, nextcls, 
	parm_no__, ehb, accp, ind)
integer *startx, *starty, *startz, *lastx, *lasty, *lastz;
shortint *space_grid__;
integer *clshd, *clsatm;
logical *excluded;
integer *search_mode__, *donp, *nextcls;
shortint *parm_no__;
doublereal *ehb;
integer *accp, *ind;
{
    /* System generated locals */
    integer space_grid_dim1, space_grid_dim2, space_grid_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer iatm;
    extern /* Subroutine */ int hbenergy_();
    static integer jx, jy, jz, ispace, ihbond, jhbond, khbond, psn;

/* #include "params.inc" */
/* #include "grid.inc" */
/* #include "cg.inc" */
/* #include "coords.inc" */
/* #include "engpar.inc" */
/* #include "pstruct.inc" */
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

/*     Data structure for space grid. Pointer variables are prefixed */
/*     by GRID_ in order to avoid naming conflicts. */


/* ! */
/* !     Side chain construction options */
/* ! */
/*    Removed FORTRAN-90 style data in variable declarations ACRM 07.01.94
*/
/*     Also removed BYTE definitions (-> CHARACTER) */
/* 07.01.94      INTEGER FIRST,INDEPENDENT,ALL,COMBINATION,ITERATIVE */
/* ! */
/* !     Evaluation options */
/* ! */
/* 07.01.94      INTEGER EO_ENERGY,EO_RMS */
/* ! */
/* !     evaluate_f options */
/* ! */
/* 07.01.94      INTEGER EVAL_CODE_RMS,EVAL_CODE_ENERGY,EVAL_CODE_USER */
/* ! */
/* 07.01.94      INTEGER ST_SIDEOPT(SIDEOPTWIDTH,NSIDEOPT), */
/* 07.01.94     2        L_SIDEOPT(NSIDEOPT) */

/* OMUPD drm 20/11/91 declare as integers */

/* 07.01.94      INTEGER ST_EVALOPT(EVALOPTWIDTH,NEVALOPT) */
/* 07.01.94      INTEGER L_EVALOPT(NEVALOPT) */
/* ! */
/* ! */
/* Change name to avoid conflict */
/* 07.01.94      DATA FIRST,INDEPENDENT,ALL,COMBINATION,ITERATIVE */
/* 07.01.94     2     /1,2,3,4,5/ */
/* 07.01.94      DATA EO_ENERGY,EO_RMS/1,2/ */
/* 07.01.94      DATA EVAL_CODE_RMS,EVAL_CODE_ENERGY,EVAL_CODE_USER/1,2,3/
 */
/* 07.01.94      DATA ST_SIDEOPT,L_SIDEOPT */
/* 07.01.94     2     /'FIRS','T',' ', */
/* 07.01.94     3     'INDE','PEND','ENT', */
/* 07.01.94     4     'ALL',' ',' ', */
/* 07.01.94     5     'COMB','INAT','ION', */
/* 07.01.94     6     'ITER','ATIV','E', */
/* 07.01.94     7     5,11,3,11,9/ */
/* 07.01.94      DATA ST_EVALOPT,L_EVALOPT/'ENER','GY','RMS',' ',6,3/ */
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
    --accp;
    parm_no__ -= 101;
    --nextcls;
    --donp;
    --excluded;
    --clsatm;
    --clshd;
    space_grid_dim1 = grid_1.ngridx;
    space_grid_dim2 = grid_1.ngridy;
    space_grid_offset = space_grid_dim1 * (space_grid_dim2 + 1) + 1;
    space_grid__ -= space_grid_offset;

    /* Function Body */
    i__1 = *lastx;
    for (jx = *startx; jx <= i__1; ++jx) {
	i__2 = *lasty;
	for (jy = *starty; jy <= i__2; ++jy) {
	    i__3 = *lastz;
	    for (jz = *startz; jz <= i__3; ++jz) {
		ispace = space_grid__[jx + (jy + jz * space_grid_dim2) * 
			space_grid_dim1];
		if (ispace != 0) {
		    psn = clshd[ispace];
L5:
		    if (psn != 0) {
			iatm = clsatm[psn];
			if (! excluded[iatm]) {
			    if (*search_mode__ == 3 && accp[iatm] != 0) {
				ihbond = pstruct_1.hbdonr[donp[*ind] - 1];
				jhbond = iatm;
				khbond = *ind;
/*                             EVALUATE-ENERGY
-FOR-ONE-HYDROGEN-BOND */
				hbenergy_(&ihbond, &jhbond, &khbond, &
					parm_no__[101], ehb);
			    } else if (*search_mode__ == 4 && donp[iatm] != 0)
				     {
				ihbond = pstruct_1.hbdonr[donp[iatm] - 1];
				jhbond = *ind;
				khbond = iatm;
/*                             EVALUATE-ENERGY
-FOR-ONE-HYDROGEN-BOND */
				hbenergy_(&ihbond, &jhbond, &khbond, &
					parm_no__[101], ehb);
			    }
			}
			psn = nextcls[psn];
			goto L5;
		    }
		}
/* L20: */
	    }
/* L50: */
	}
/* L100: */
    }
    return 0;
} /* rnhbndsrch_ */

