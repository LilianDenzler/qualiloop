/* srchnat1.f -- translated by f2c (version of 23 April 1993  18:34:30).
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
    integer dbg_cgen__, dbg_clschn__, dbg_alloc__, dbg_allstk__, dbg_allhp__;
} dbg_;

#define dbg_1 dbg_

struct {
    shortint hbdon[4100], hbacc[4100];
    real hbcut, hbacut;
} hbonds_;

#define hbonds_1 hbonds_

struct {
    integer natoms, nres, nsegs, nbonds, nangs, nptors, nitors, nhbs, nnbs, 
	    ndonat, naccat, nbpar, napar, nptpar, nitpar, nhbpar, natyps, 
	    nbauto;
} values_;

#define values_1 values_

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
    shortint bndat1[4000]	/* was [100][40] */, bndat2[4000]	/* 
	    was [100][40] */, angat1[6000]	/* was [150][40] */, angat2[
	    6000]	/* was [150][40] */, angat3[6000]	/* was [150][
	    40] */, torat1[1400]	/* was [35][40] */, torat2[1400]	
	    /* was [35][40] */, torat3[1400]	/* was [35][40] */, torat4[
	    1400]	/* was [35][40] */, impat1[880]	/* was [22][40] */, 
	    impat2[880]	/* was [22][40] */, impat3[880]	/* was [22][40] */, 
	    impat4[880]	/* was [22][40] */, exclnb[5600]	/* was [140][
	    40] */, nexcld[2800]	/* was [70][40] */, acpthb[320]	/* 
	    was [8][40] */, aan1hb[320]	/* was [8][40] */, aan2hb[320]	/* 
	    was [8][40] */, donrhb[320]	/* was [8][40] */, dan1hb[320]	/* 
	    was [8][40] */, dan2hb[320]	/* was [8][40] */, dhydhb[320]	/* 
	    was [8][40] */, fstgrp[400]	/* was [10][40] */, lstgrp[400]	/* 
	    was [10][40] */, acindx[2800]	/* was [70][40] */, nparam[
	    400]	/* was [10][40] */, nparus[400]	/* was [10][40] */, 
	    nparmx[10];
    integer nreses;
    real atmmas[100], atmchg[70], b1bld[2800]	/* was [70][40] */, b2bld[
	    2800]	/* was [70][40] */, a1bld[2800]	/* was [70][40] */, 
	    a2bld[2800]	/* was [70][40] */, torbld[2800]	/* was [70][
	    40] */, bldat1[2800]	/* was [70][40] */, bldat2[2800]	
	    /* was [70][40] */, bldat3[2800]	/* was [70][40] */, bldat4[
	    2800]	/* was [70][40] */, namres[40], acodes[100], namatm[
	    2800]	/* was [70][40] */, grpnam[400]	/* was [10][40] */;
} restop_;

#define restop_1 restop_

/* Subroutine */ int srchnat1_(ind, space_grid__, search_mode__, maxevdw, 
	ignore_evdw__, cntnbx, nbxa, excluded, parm_no__, clshd, nextcls, 
	clsatm, donp, accp, qside, sidehits, resbya, impact, eel, evdw, ehb, 
	outofbound)
integer *ind;
shortint *space_grid__;
integer *search_mode__;
real *maxevdw;
logical *ignore_evdw__;
integer *cntnbx;
shortint *nbxa;
logical *excluded;
shortint *parm_no__;
integer *clshd, *nextcls, *clsatm, *donp, *accp;
logical *qside;
integer *sidehits, *resbya;
logical *impact;
doublereal *eel, *evdw, *ehb;
logical *outofbound;
{
    /* System generated locals */
    integer space_grid_dim1, space_grid_dim2, space_grid_offset, nbxa_dim1, 
	    nbxa_offset, i__1;
    real r__1;

    /* Local variables */
    static integer irad;
    static real xind, yind, zind;
    extern /* Subroutine */ int rnnbsrch_();
    static integer i, lastx, lasty, lastz, ix, iy, iz, startx, starty, startz;
    extern integer chmceil_();
    extern /* Subroutine */ int compspcgrd_(), rnhbndsrch_();
    static real cut2;


/*     Performs a search of the space near atom, IND. The actual */
/*     operation is determined by the variable, SEARCH_MODE, which is */
/*     decoded as follows: */

/*     SEARCH_MODE value       Function */

/*     CONTACT                 Search for contacts which exceed MAXEVDW */
/*                             in energy */
/*     NBOND_ENERGY            Evaluate non-bonded energy terms. */
/*     DONOR_HBOND_ENERGY      Calculate hydrogen bond energy for */
/*                             hydrogen atom, IND. */
/*     ACCEPTOR_HBOND_ENERGY   Calculate hydrogen bond energy for the */
/*                             hydrogen bond acceptor, IND. */

/*     If the atom, IND, does not fit in the grid, OUTOFBOUND will */
/*     be returned as true so that the grid can be enlarged. */

/* #include "grid.inc" */
/* #include "params.inc" */
/* #include "cg.inc" */
/* #include "coords.inc" */
/* #include "dbg.inc" */
/* #include "hbonds.inc" */
/* #include "values.inc" */
/* #include "engpar.inc" */
/* #include "pstruct.inc" */
/* #include "restop.inc" */

/*     Data structure for space grid. Pointer variables are prefixed */
/*     by GRID_ in order to avoid naming conflicts. */


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
/* *      NAME: HBONDS                                                   *
 */
/* *  FUNCTION: To declare hydrogen bonding parameters                   *
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
/*  Variable   Array Bounds              Description */
/*  --------   ------------              ----------- */

/*  HBCUT          -          The hydrogen bond distance limit */
/*  HBACUT         -          The hydrogen bond angle limit */
/*  HBDON        MAXHB        The heavy atom donor for each hydrogen bond 
*/
/*  HBACC        MAXHB        The heavy atom acceptor for each hydrogen */

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
/* Argument declarations */
/* ***********************************************************************
 */
/* *      NAME: RESTOP                                                   *
 */
/* *  FUNCTION: To declare the variable arrays used to store the         *
 */
/* *            information read in from the main residue topology file. *
 */
/* *            These are then used to derive the potential energy       *
 */
/* *            function parameters                                      *
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
/* Variable   Array Bounds                 Description */
/* --------   ------------                 ----------- */
/* Parameter variables: */
/* MXATMR         -        The maximum number of atoms per residue */
/* MXBNDR         -        The maximum number of bonds per residue */
/* MXANGR         -        The maximum number of bond angles per residue 
*/
/* MXTORR         -        The maximum number of torsions per residue */
/* MXIMPR         -        The maximum number of improper torsions per */
/*                         residue */
/* MXNBER         -        The maximum number of non-bond exclusions */
/*                         per residue */
/* MXHBAR         -        The maximum number of H-bond acceptors per */
/*                         residue */
/* MXHBDR         -        The maximum number of H-bond donors per */
/*                         residue */
/* MXICR          -        The maximum number of internal coordinates */
/*                         per residue */
/* MXGRPR         -        The maximum number of groups per residue */
/* MAXIMR         -        The maximum allowed number of residues */
/* NINFO          -        Total number of different parameter types */
/* Atom information variables: */
/* ATMMAS      (MAXATT)    List of atomic masses, indexed with ACODES */
/* ATMCHG  (MXATMR,MAXIMR) List of atomic charges */
/* ACODES      (MAXATT)    Atom type names used with ARMASS */
/* ACINDX  (MXATMR,MAXIMR) Index of atom type codes for ACODES */
/* NAMATM  (MXATMR,MAXIMR) List of IUPAC atom names */
/* Group information variables: */
/* FSTGRP  (MXGRPR,MAXIMR) List of first atoms in each group */
/* GRPNAM  (MXGRPR,MAXIMR) List of group names for each residue */
/* LSTGRP  (MXGRPR,MAXIMR) List of last atoms in each group */
/* Residue information variables: */
/* NRESES         -        The number of residues in use */
/* NAMRES      (MAXIMR)    List of residue names */
/* Connectivity information: */
/* BNDAT1  (MXBNDR,MAXIMR) Atom number 1 in each bond */
/* BNDAT2  (MXBNDR,MAXIMR) Atom number 2 in each bond */
/* ANGAT1  (MXANGR,MAXIMR) Atom number 1 in each bond angle */
/* ANGAT2  (MXANGR,MAXIMR) Atom number 2 in each bond angle */
/* ANGAT3  (MXANGR,MAXIMR) Atom number 3 in each bond angle */
/* TORAT1  (MXTORR,MAXIMR) Atom number 1 in each proper torsion */
/* TORAT2  (MXTORR,MAXIMR) Atom number 2 in each proper torsion */
/* TORAT3  (MXTORR,MAXIMR) Atom number 3 in each proper torsion */
/* TORAT4  (MXTORR,MAXIMR) Atom number 4 in each proper torsion */
/* IMPAT1  (MXIMPR,MAXIMR) Atom number 1 in each improper torsion */
/* IMPAT2  (MXIMPR,MAXIMR) Atom number 2 in each improper torsion */
/* IMPAT3  (MXIMPR,MAXIMR) Atom number 3 in each improper torsion */
/* IMPAT4  (MXIMPR,MAXIMR) Atom number 4 in each improper torsion */
/* EXCLNB  (MXNBER,MAXIMR) List non-bonded exclusions. This is indexed */
/*                         for each atom with NEXCLD */
/* NEXCLD  (MXATMR,MAXIMR) The number of non-bond exclusions per atom */
/* ACPTHB  (MXHBAR,MAXIMR) List of hydrogen bond acceptor atoms */
/* AAN1HB  (MXHBAR,MAXIMR) First antecedent to ACPTHB */
/* AAN2HB  (MXHBAR,MAXIMR) Second antecedent to ACPTHB */
/* DONRHB  (MXHBDR,MAXIMR) List of hydrogen bond donor atoms */
/* DAN1HB  (MXHBDR,MAXIMR) First antecedent to DONRHB */
/* DAN2HB  (MXHBDR,MAXIMR) Second antecedent to DONRHB */
/* DHYDHB  (MXHBDR,MAXIMR) Hydrogens attached to DONRHB in H-bond */
/* NPARAM  (NINFO,MAXIMR)  The total number of some system parameters. */
/*                         atoms, bonds, angles, torsions, improper */
/*                         torsions, non-bond exclusions, H-bond donors, 
*/
/*                         H-bond acceptors, internal coordinates and */
/*                         groups */
/* NPARUS      (NINFO)     As NPARAM but this gives the number which are 
*/
/*                         actually to be used. */
/* NPARMX      (NINFO)     Array bounds for NPARAM */
/* Build information: */
/*BLDAT1  (MXICR,MAXIMR)  First atom for each internal coordinate in build
*/
/*BLDAT2  (MXICR,MAXIMR)  Second atom for each internal coordinate in buil
d*/
/*BLDAT3  (MXICR,MAXIMR)  Third atom for each internal coordinate in build
*/
/*BLDAT4  (MXICR,MAXIMR)  Fourth atom for each internal coordinate in buil
d*/
/* B1BLD   (MXICR,MAXIMR)  First bond length for internal coordinates */
/* B2BLD   (MXICR,MAXIMR)  Next bond length for internal coordinates */
/* A1BLD   (MXICR,MAXIMR)  First bond angle for internal coordinates */
/* A2BLD   (MXICR,MAXIMR)  Next bond angle for internal coordinates */
/* TORBLD  (MXICR,MAXIMR)  Torsion angle value for internal coordinates */

/* Parameter declarations */
/* Other declarations */
/* Local declarations */

    /* Parameter adjustments */
    --resbya;
    --sidehits;
    --qside;
    --accp;
    --donp;
    --clsatm;
    --nextcls;
    --clshd;
    parm_no__ -= 101;
    --excluded;
    nbxa_dim1 = grid_1.maxnbx;
    nbxa_offset = nbxa_dim1 + 1;
    nbxa -= nbxa_offset;
    --cntnbx;
    space_grid_dim1 = grid_1.ngridx;
    space_grid_dim2 = grid_1.ngridy;
    space_grid_offset = space_grid_dim1 * (space_grid_dim2 + 1) + 1;
    space_grid__ -= space_grid_offset;

    /* Function Body */
    if (*search_mode__ == 1) {
	*impact = FALSE_;
    }
    compspcgrd_(ind, &ix, &iy, &iz, outofbound);
    if (*outofbound) {
	return 0;
    }
    if (*search_mode__ == 1) {
	irad = 1;
	cut2 = (float)1.;
    } else if (*search_mode__ == 2) {
/* ACRM */
/*    NBCUT taken of the heap and simply placed in */
/*    the NBND common block */
/* CC            IRAD = CHMCEIL(RHEAP(BNBND(NBCUT))*RECIPGRID) */
/* CC            CUT2 = (RHEAP(BNBND(NBCUT))*RECIPGRID)**2 */
	r__1 = engpar_1.nbcut * grid_1.recipgrid;
	irad = chmceil_(&r__1);
/* Computing 2nd power */
	r__1 = engpar_1.nbcut * grid_1.recipgrid;
	cut2 = r__1 * r__1;
    } else if (*search_mode__ == 4) {
	r__1 = hbonds_1.hbcut * grid_1.recipgrid;
	irad = chmceil_(&r__1);
/* Computing 2nd power */
	r__1 = hbonds_1.hbcut * grid_1.recipgrid;
	cut2 = r__1 * r__1;
    } else if (*search_mode__ == 3) {
	r__1 = hbonds_1.hbcut * grid_1.recipgrid;
	irad = chmceil_(&r__1);
/* Computing 2nd power */
	r__1 = hbonds_1.hbcut * grid_1.recipgrid;
	cut2 = r__1 * r__1;
    }
    startx = ix - irad;
    starty = iy - irad;
    startz = iz - irad;
    if (startx < 1) {
	startx = 1;
    }
    if (starty < 1) {
	starty = 1;
    }
    if (startz < 1) {
	startz = 1;
    }
    lastx = ix + irad;
    lasty = iy + irad;
    lastz = iz + irad;
    if (lastx > grid_1.ngridx) {
	lastx = grid_1.ngridx;
    }
    if (lasty > grid_1.ngridy) {
	lasty = grid_1.ngridy;
    }
    if (lastz > grid_1.ngridz) {
	lastz = grid_1.ngridz;
    }
    i__1 = cntnbx[*ind];
    for (i = 1; i <= i__1; ++i) {
	excluded[nbxa[i + *ind * nbxa_dim1]] = TRUE_;
/* L100: */
    }
    xind = coords_1.xcart[*ind - 1];
    yind = coords_1.ycart[*ind - 1];
    zind = coords_1.zcart[*ind - 1];
    if (*search_mode__ == 2 || *search_mode__ == 1) {
	rnnbsrch_(search_mode__, &startx, &lastx, &starty, &lasty, &startz, &
		lastz, &space_grid__[space_grid_offset], &clshd[1], &clsatm[1]
		, &excluded[1], &parm_no__[101], maxevdw, impact, ind, &
		resbya[1], &qside[1], &sidehits[1], eel, ignore_evdw__, evdw, 
		&nextcls[1], &xind, &yind, &zind);
    } else {
	rnhbndsrch_(&startx, &starty, &startz, &lastx, &lasty, &lastz, &
		space_grid__[space_grid_offset], &clshd[1], &clsatm[1], &
		excluded[1], search_mode__, &donp[1], &nextcls[1], &parm_no__[
		101], ehb, &accp[1], ind);
    }
    i__1 = cntnbx[*ind];
    for (i = 1; i <= i__1; ++i) {
	excluded[nbxa[i + *ind * nbxa_dim1]] = FALSE_;
/* L200: */
    }
    return 0;
} /* srchnat1_ */

