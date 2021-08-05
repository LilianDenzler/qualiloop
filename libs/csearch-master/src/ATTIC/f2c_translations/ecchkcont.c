/* ecchkcont.f -- translated by f2c (version of 23 April 1993  18:34:30).
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
    real xcart[6150], ycart[6150], zcart[6150], xwork[6150], ywork[6150], 
	    zwork[6150];
} coords_;

#define coords_1 coords_

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

static integer c__1 = 1;
static integer c__3 = 3;

/* Subroutine */ int ecchkcont_(iatm, atomno, ori, ecx, ecy, ecz, avoidx, 
	avoidr, avoidphi, rcontact, startphi, lastphi, impact, sidehits, 
	resbya, cntnbx, nbxa, qside)
integer *iatm, *atomno;
real *ori, *ecx, *ecy, *ecz, *avoidx, *avoidr, *avoidphi, *rcontact, *
	startphi, *lastphi;
logical *impact;
integer *sidehits, *resbya, *cntnbx;
shortint *nbxa;
logical *qside;
{
    /* Format strings */
    static char fmt_9001[] = "(\002 EC_CHECK_CONTACT: \002,a20,\002 and \002\
,a20)";
    static char fmt_9002[] = "(\002           Angles: \002,2f10.3,\002 degre\
es.\002)";
    static char fmt_9003[] = "(\002 Origin is at \002,3f8.3,\002. Clump atom\
: \002,3f8.3,\002.\002)";

    /* System generated locals */
    integer nbxa_dim1, nbxa_offset, i__1;
    real r__1, r__2, r__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(), acos(), atan2();
    integer s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    static real c, d;
    static integer i, i1;
    static real delta1;
    extern /* Subroutine */ int forprtatm_();
    static real aa, bb, cc, dd, deltap;
    static char buffer[100];
    static real xs, xt, yt, zt, ys, zs;
    static char st1[21], st2[21];
    extern /* Subroutine */ int cprint_();

    /* Fortran I/O blocks */
    static icilist io___20 = { 0, buffer, 0, fmt_9001, 100, 1 };
    static icilist io___21 = { 0, buffer, 0, fmt_9002, 100, 1 };
    static icilist io___22 = { 0, buffer, 0, fmt_9003, 100, 1 };



/*     Determines the range of torsion angles for the free torsion in the 
*/
/*     clump of which atom number, ATOMNO, is a member, where there is no 
*/
/*     contact closer than those found in RCONTACT to the atom numbered */
/*     IATM. ORI, ECX, ECY, and ECZ define the local coordinate system */
/*     for the clump. AVOIDX gives the distance along the x axis for */
/*     ATOMNO, AVOIDR gives the rotation radius for that atom, AVOIDPHI */
/*     gives the offset from the free torsion, RCONTACT is an array of */
/*     van der Waals radii for all atom types with ATOMNO. STARTPHI and */
/*     LASTPHI give the allowed range where STARTPHI is always less than 
*/
/*     or equal to LASTPHI and STARTPHI is within [-pi,pi]. If no */
/*     placement is possible, then STARTPHI is returned with -10.0. If */
/*     all placements are possible, then IMPACT is turned off. The reason 
*/
/*     for this redundancy over STARTPHI and LASTPHI is the imprecision */
/*     of floating point numbers. SIDEHITS records any sidechain hits for 
*/
/*     ATOMNO in the event that both IATM and ATOMNO are sidechain atoms 
*/
/*     and IATM could make close contact. Finally, RESBYA gives residue */
/*     numbers by atom to permit quick access into SIDEHITS. */

/*      IMPLICIT INTEGER(A-Z) */
/* OM OMUPD BNJ 2/9/91 */
/* OM */
/*     Moved variable definitions before INCLUDES (since these define */
/*     data) ACRM 07.01.94 */

/* OM   ACRM Change BYTE to CHARACTER */
/* OM      BYTE ST1(21),ST2(21) */
/* #include "params.inc" */
/* #include "cg.inc" */
/* #include "values.inc" */
/* #include "coords.inc" */
/* #include "dbg.inc" */
/* #include "grid.inc" */
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

/*     Data structure for space grid. Pointer variables are prefixed */
/*     by GRID_ in order to avoid naming conflicts. */


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
    --qside;
    nbxa_dim1 = grid_1.maxnbx;
    nbxa_offset = nbxa_dim1 + 1;
    nbxa -= nbxa_offset;
    --cntnbx;
    --resbya;
    --sidehits;
    --rcontact;
    --ecz;
    --ecy;
    --ecx;
    --ori;

    /* Function Body */
    *startphi = (float)-3.141592653589794;
    *lastphi = (float)3.141592653589794;
    *impact = FALSE_;
    i__1 = cntnbx[*atomno];
    for (i = 1; i <= i__1; ++i) {
	if (nbxa[i + *atomno * nbxa_dim1] == *iatm) {
	    return 0;
	}
/* L100: */
    }

/*     Get coordinates of IATM in coordinate system for ATOMNO. */

    xt = coords_1.xcart[*iatm - 1] - ori[1];
    yt = coords_1.ycart[*iatm - 1] - ori[2];
    zt = coords_1.zcart[*iatm - 1] - ori[3];
    xs = xt * ecx[1] + yt * ecx[2] + zt * ecx[3];
    ys = xt * ecy[1] + yt * ecy[2] + zt * ecy[3];
    zs = xt * ecz[1] + yt * ecz[2] + zt * ecz[3];
    i1 = engpar_1.atflag[pstruct_1.atcode[*iatm - 1] - 1];
/* Computing 2nd power */
    r__1 = rcontact[i1];
/* Computing 2nd power */
    r__2 = *avoidx - xs;
    aa = r__1 * r__1 - r__2 * r__2;
/* Computing 2nd power */
    r__1 = *avoidr;
/* Computing 2nd power */
    r__2 = ys;
/* Computing 2nd power */
    r__3 = zs;
    bb = aa - r__1 * r__1 - r__2 * r__2 - r__3 * r__3;
    c = -(doublereal)bb / (float)2. / *avoidr;
    cc = c * c;
/* Computing 2nd power */
    r__1 = ys;
/* Computing 2nd power */
    r__2 = zs;
    dd = r__1 * r__1 + r__2 * r__2;
    if (dd < cc && c > (float)0.) {
    } else if (dd <= cc && c <= (float)0.) {
	*startphi = (float)-10.;
	*impact = TRUE_;
/*        CHECK-SIDEHITS */
	if (qside[*iatm] && qside[*atomno]) {
	    if (dbg_1.dbg_cgen__ > 1) {
		cprint_(" Next contact entered in sidehits array.", 40L);
	    }
	    sidehits[resbya[*iatm]] = *atomno;
	}
    } else {
	*impact = TRUE_;
	d = sqrt(dd);
	delta1 = acos(c / d);
	deltap = atan2(zs, ys) - *avoidphi;
	*startphi = deltap + delta1;
	*lastphi = deltap - delta1 + 6.283185307179588;
L150:
	if (*startphi > 3.141592653589794) {
	    *startphi += -6.283185307179588;
	    *lastphi += -6.283185307179588;
	    goto L150;
	}
L200:
	if (*startphi < -3.141592653589794) {
	    *startphi += 6.283185307179588;
	    *lastphi += 6.283185307179588;
	    goto L200;
	}
/*        CHECK-SIDEHITS */
	if (qside[*iatm] && qside[*atomno]) {
	    if (dbg_1.dbg_cgen__ > 1) {
		cprint_(" Next contact entered in sidehits array.", 40L);
	    }
	    sidehits[resbya[*iatm]] = *atomno;
	}
    }
    if (dbg_1.dbg_cgen__ > 1) {
	if (dbg_1.dbg_cgen__ > 2 || *startphi != -3.141592653589794 || *
		lastphi != 3.141592653589794) {
/* OM ACRM Changed to use FORPRTATM rather then BYTES... */
	    forprtatm_(st1, iatm, 21L);
	    forprtatm_(st2, atomno, 21L);
	    s_wsfi(&io___20);
	    do_fio(&c__1, st1, 21L);
	    do_fio(&c__1, st2, 21L);
	    e_wsfi();
	    cprint_(buffer, 100L);
	    s_wsfi(&io___21);
	    d__1 = *startphi / .017453292519943299;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    d__2 = *lastphi / .017453292519943299;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    cprint_(buffer, 100L);
	}
	if (dbg_1.dbg_cgen__ > 2) {
	    s_wsfi(&io___22);
	    do_fio(&c__3, (char *)&ori[1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&coords_1.xcart[*iatm - 1], (ftnlen)sizeof(
		    real));
	    do_fio(&c__1, (char *)&coords_1.ycart[*iatm - 1], (ftnlen)sizeof(
		    real));
	    do_fio(&c__1, (char *)&coords_1.zcart[*iatm - 1], (ftnlen)sizeof(
		    real));
	    e_wsfi();
	    cprint_(buffer, 100L);
	}
    }
    return 0;
} /* ecchkcont_ */

