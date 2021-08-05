/* rnnbsrch.f -- translated by f2c (version of 23 April 1993  18:34:30).
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

struct {
    integer natoms, nres, nsegs, nbonds, nangs, nptors, nitors, nhbs, nnbs, 
	    ndonat, naccat, nbpar, napar, nptpar, nitpar, nhbpar, natyps, 
	    nbauto;
} values_;

#define values_1 values_

/* Table of constant values */

static integer c__1 = 1;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int rnnbsrch_(search_mode__, startx, lastx, starty, lasty, 
	startz, lastz, space_grid__, clshd, clsatm, excluded, parm_no__, 
	maxevdw, impact, ind, resbya, qside, sidehits, eel, ignore_evdw__, 
	evdw, nextcls, xind, yind, zind)
integer *search_mode__, *startx, *lastx, *starty, *lasty, *startz, *lastz;
shortint *space_grid__;
integer *clshd, *clsatm;
logical *excluded;
shortint *parm_no__;
real *maxevdw;
logical *impact;
integer *ind, *resbya;
logical *qside;
integer *sidehits;
doublereal *eel;
logical *ignore_evdw__;
doublereal *evdw;
integer *nextcls;
real *xind, *yind, *zind;
{
    /* Format strings */
    static char fmt_9001[] = "(\002 Close contact failed. R =\002,f5.2,\002 \
E = \002,1pg12.5)";
    static char fmt_9002[] = "(\002 Atoms:\002,4(1x,a4),\002 and\002,4(1x,a4\
))";
    static char fmt_9003[] = "(\002 Sidechain collision, atoms:\002,2(1x,a4,\
i5))";
    static char fmt_9004[] = "(\002 Energy for atoms:\002,2(1x,a4,i5),\002 R\
 =\002,f5.2,\002 EVDW = \002,1pg12.5,\002 ELEC = \002,1pg12.5)";

    /* System generated locals */
    integer space_grid_dim1, space_grid_dim2, space_grid_offset, i__1, i__2, 
	    i__3;
    real r__1;

    /* Builtin functions */
    integer s_wsfi();
    double sqrt();
    integer do_fio(), e_wsfi();

    /* Local variables */
    static integer iatm;
    static real delx, dely, delz, cgind, r2, r4, r6, f1, f2;
    static integer ic;
    static real r12;
    static integer itcind, offind, jx;
    static real cut2_search__;
    static integer jy, jz, ispace;
    static real enb, ss;
    static integer indseg;
    extern integer getseg_();
    static char buffer[100];
    extern /* Subroutine */ int cprint_();
    static integer psn, iatmseg, itciatm;
    static real eel1;

    /* Fortran I/O blocks */
    static icilist io___27 = { 0, buffer, 0, fmt_9001, 100, 1 };
    static icilist io___28 = { 0, buffer, 0, fmt_9002, 100, 1 };
    static icilist io___29 = { 0, buffer, 0, fmt_9003, 100, 1 };
    static icilist io___31 = { 0, buffer, 0, fmt_9004, 100, 1 };


/* #include "params.inc" */
/* #include "grid.inc" */
/* #include "cg.inc" */
/* #include "coords.inc" */
/* #include "engpar.inc" */
/* #include "pstruct.inc" */
/* #include "dbg.inc" */
/* #include "values.inc" */
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
    --nextcls;
    --sidehits;
    --qside;
    --resbya;
    parm_no__ -= 101;
    --excluded;
    --clsatm;
    --clshd;
    space_grid_dim1 = grid_1.ngridx;
    space_grid_dim2 = grid_1.ngridy;
    space_grid_offset = space_grid_dim1 * (space_grid_dim2 + 1) + 1;
    space_grid__ -= space_grid_offset;

    /* Function Body */
    if (*search_mode__ == 1) {
	cut2_search__ = cg_1.grid2;
    } else {
	cut2_search__ = cg_1.cutnb2;
    }
    cgind = pstruct_1.atchrg[*ind - 1];
    itcind = engpar_1.atflag[pstruct_1.atcode[*ind - 1] - 1];
    offind = cg_1.ioff[itcind - 1];
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
			    delx = coords_1.xcart[iatm - 1] - *xind;
			    dely = coords_1.ycart[iatm - 1] - *yind;
			    delz = coords_1.zcart[iatm - 1] - *zind;
			    ss = delx * delx + dely * dely + delz * delz;
			    if (ss < cut2_search__) {
				if (*search_mode__ == 1) {
				    r2 = (float)1. / ss;
				    r4 = r2 * r2;
				    r6 = r4 * r2;
				    r12 = r6 * r6;
				    itciatm = engpar_1.atflag[
					    pstruct_1.atcode[iatm - 1] - 1];
				    ic = parm_no__[itcind + itciatm * 100];
				    f1 = engpar_1.vdwr12[ic - 1] * r12;
				    f2 = engpar_1.vdwr6[ic - 1] * r6;
				    enb = f1 - f2;
				    if (enb > *maxevdw) {
					*impact = TRUE_;
					if (dbg_1.dbg_cgen__ > 1) {
					    indseg = getseg_(&resbya[*ind], 
						    pstruct_1.segndx, &
						    values_1.nsegs);
					    iatmseg = getseg_(&resbya[iatm], 
						    pstruct_1.segndx, &
						    values_1.nsegs);
					    s_wsfi(&io___27);
					    r__1 = sqrt(ss);
					    do_fio(&c__1, (char *)&r__1, (
						    ftnlen)sizeof(real));
					    do_fio(&c__1, (char *)&enb, (
						    ftnlen)sizeof(real));
					    e_wsfi();
					    cprint_(buffer, 100L);
					    s_wsfi(&io___28);
					    do_fio(&c__1, (char *)&
						    pstruct_1.segid[indseg - 
						    1], (ftnlen)sizeof(real));
					    do_fio(&c__1, (char *)&
						    pstruct_1.resnme[resbya[*
						    ind] - 1], (ftnlen)sizeof(
						    real));
					    do_fio(&c__1, (char *)&
						    pstruct_1.resid[resbya[*
						    ind] - 1], (ftnlen)sizeof(
						    real));
					    do_fio(&c__1, (char *)&
						    pstruct_1.atmnme[*ind - 1]
						    , (ftnlen)sizeof(real));
					    do_fio(&c__1, (char *)&
						    pstruct_1.segid[iatmseg - 
						    1], (ftnlen)sizeof(real));
					    do_fio(&c__1, (char *)&
						    pstruct_1.resnme[resbya[
						    iatm] - 1], (ftnlen)
						    sizeof(real));
					    do_fio(&c__1, (char *)&
						    pstruct_1.resid[resbya[
						    iatm] - 1], (ftnlen)
						    sizeof(real));
					    do_fio(&c__1, (char *)&
						    pstruct_1.atmnme[iatm - 1]
						    , (ftnlen)sizeof(real));
					    e_wsfi();
					    cprint_(buffer, 100L);
					}
					if (qside[iatm] && qside[*ind]) {
					    if (dbg_1.dbg_cgen__ > 1) {
			  s_wsfi(&io___29);
			  do_fio(&c__1, (char *)&pstruct_1.atmnme[*ind - 1], (
				  ftnlen)sizeof(real));
			  do_fio(&c__1, (char *)&(*ind), (ftnlen)sizeof(
				  integer));
			  do_fio(&c__1, (char *)&pstruct_1.atmnme[iatm - 1], (
				  ftnlen)sizeof(real));
			  do_fio(&c__1, (char *)&iatm, (ftnlen)sizeof(integer)
				  );
			  e_wsfi();
			  cprint_(buffer, 100L);
					    }
					    sidehits[resbya[iatm]] = *ind;
					}
					goto L200;
				    }
				} else {
				    r2 = (float)1. / ss;
				    if (cg_1.cons_die__) {
					eel1 = cgind * (float)332.0716 * 
						pstruct_1.atchrg[iatm - 1] * 
						sqrt(r2) / cg_1.epsilon;
				    } else {
					eel1 = cgind * (float)332.0716 * 
						pstruct_1.atchrg[iatm - 1] * 
						r2;
				    }
				    *eel += eel1;
				    if (! (*ignore_evdw__)) {
					r4 = r2 * r2;
					r6 = r4 * r2;
					r12 = r6 * r6;
					itciatm = engpar_1.atflag[
						pstruct_1.atcode[iatm - 1] - 
						1];
					ic = parm_no__[itcind + itciatm * 100]
						;
					f1 = engpar_1.vdwr12[ic - 1] * r12;
					f2 = engpar_1.vdwr6[ic - 1] * r6;
					enb = f1 - f2;
					*evdw += enb;
				    }
				    if (dbg_1.dbg_cgen__ > 3) {
					s_wsfi(&io___31);
					do_fio(&c__1, (char *)&
						pstruct_1.atmnme[*ind - 1], (
						ftnlen)sizeof(real));
					do_fio(&c__1, (char *)&(*ind), (
						ftnlen)sizeof(integer));
					do_fio(&c__1, (char *)&
						pstruct_1.atmnme[iatm - 1], (
						ftnlen)sizeof(real));
					do_fio(&c__1, (char *)&iatm, (ftnlen)
						sizeof(integer));
					r__1 = sqrt(ss);
					do_fio(&c__1, (char *)&r__1, (ftnlen)
						sizeof(real));
					do_fio(&c__1, (char *)&enb, (ftnlen)
						sizeof(real));
					do_fio(&c__1, (char *)&eel1, (ftnlen)
						sizeof(real));
					e_wsfi();
					cprint_(buffer, 100L);
				    }
				}
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
L200:
    return 0;
} /* rnnbsrch_ */

