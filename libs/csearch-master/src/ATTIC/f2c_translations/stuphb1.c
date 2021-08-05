/* stuphb1.f -- translated by f2c (version of 23 April 1993  18:34:30).
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
    shortint hbdon[4100], hbacc[4100];
    real hbcut, hbacut;
} hbonds_;

#define hbonds_1 hbonds_

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

/* Table of constant values */

static integer c__0 = 0;

/* Subroutine */ int stuphb1_(donp, accp)
integer *donp, *accp;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos();

    /* Local variables */
    extern /* Subroutine */ int fill4_();
    static integer i;
    extern /* Subroutine */ int cprint_(), die_();


/*     Sets up the hydrogen bond data for CGEN. DONP and ACCP are */
/*     arrays which mark each atom as being a donor or acceptor, */
/*     respectively. Non-zero elements serve as markers. Other variables 
*/
/*     pertaining to hydorgen bonds are also set. */

/* #include "params.inc" */
/* #include "values.inc" */
/* #include "pstruct.inc" */
/* #include "hbonds.inc" */
/* #include "cg.inc" */
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

    /* Parameter adjustments */
    --accp;
    --donp;

    /* Function Body */
    fill4_(&donp[1], &values_1.natoms, &c__0);
    fill4_(&accp[1], &values_1.natoms, &c__0);
/* ACRM Do we actually need this?? */
    i__1 = values_1.ndonat;
    for (i = 1; i <= i__1; ++i) {
	if (pstruct_1.hbdhyd[i - 1] == 0) {
	    cprint_("Error in CGEN -- Hydrogen bonding hydrogens must be pre\
sent.", 60L);
	    die_();
	}
	donp[pstruct_1.hbdhyd[i - 1]] = i;
/* L100: */
    }
    i__1 = values_1.naccat;
    for (i = 1; i <= i__1; ++i) {
	accp[pstruct_1.hbacpt[i - 1]] = i;
/* L200: */
    }
    cg_1.cuthb2 = hbonds_1.hbcut * hbonds_1.hbcut;
    cg_1.cos_cuthba__ = cos(hbonds_1.hbacut * .017453292519943299);
    return 0;
} /* stuphb1_ */

