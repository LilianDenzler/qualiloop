/* clschna.f -- translated by f2c (version of 23 April 1993  18:34:30).
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
    integer dbg_cgen__, dbg_clschn__, dbg_alloc__, dbg_allstk__, dbg_allhp__;
} dbg_;

#define dbg_1 dbg_

/* Table of constant values */

static integer c__1 = 1;
static integer c__9 = 9;
static integer c__3 = 3;
static real c_b352 = (float)100.;
static real c_b366 = (float)1.;
static integer c__4 = 4;

/* OM OMUPD BNJ 19/11/91 */
/* OM */
/* OM CISTRANS ==> ENTZUS */
/* OM PEPOMEGA ==> POMEGA */
/* OM NEWANGLE ==> NANGLE */
/* OM ITERMAX ==> ITMAX */
/* OM ITERPHI ==> ITRPHI */
/* OM DIFFPHI ==> DIFPHI */
/* OM GLOWPHI ==> GLOPHI */
/* OM BASEMAXIT ==> BSMXIT */
/* OM VARDEBUG ==> VARDBG */
/* OM MAXITRF ==> MXITRF */
/* OM OLD_NP ==> NPOLD */
/* OM SMALLEST_NP ==> NPLITL */
/* OM REFER_NP ==> REFNP */
/* OM NEW_NP ==> NPNEW */
/* OM SMALLEST_SGNS ==> SGNSML */
/* OM REFER_SGNS ==> SGNREF */
/* OM NPOLDHI1 ==> NPOLH1 */
/* OM SAVED_I* ==> SAVI* */
/* OM MODIFY_TRIES ==> MODTRI */
/* OM NSDSTEP ==> NSDSTP */
/* OM MAXITMIN ==> MXITMN */
/* OM MAX_MODIFY_TRIES ==> MXMODT */
/* OM BIG_G ==> BIGGEE */
/* OM SMALLEST_G ==> SMLGEE */
/* OM SMALLEST_OMEGA1 ==> SMLOM1 */
/* OM OLD_PHI1 ==> OLDP1 */
/* OM PHIRANGE ==> PHIRNG */
/* OM PHI1_FRAC ==> P1FRAC */
/* OM REFER_OMEGA1 ==> REFOM1 */
/* OM DISPLAY_SG ==> DISPSG */
/* OM GPHIPREV ==> GPPREV */
/* OM SAVED_ANGLE ==> SAVANG */
/* OM DOTMOVEGRAD ==> DTMVGD */
/* OM NORMSHIFT ==> NMSHFT */
/* OM LENSHIFT ==> LNSHFT */
/* OM PREVANGLE ==> PRVANG */
/* OM LENMOVE ==> LENMV */
/* OM POS_SEEN ==> POSSEN */
/* OM NEG_SEEN ==> NEGSEN */
/* OM ZERO_G_FOUND ==> FND0G */
/* OM GBOUND_OK ==> GBNDOK */
/* OM MODIFY_DONE ==> MODDUN */
/* OM OLDANGLE ==> OLDANG */
/* OM OLD_OMEGA1 ==> OLDOM1 */
/* OM STOP_SD ==> STOPSD */
/* OM DROP_STEP ==> DRPSTP */
/* OM GO_ON ==> CARYON */
/* OM ADJUSTED_THETA ==> ADTHET */
/* OM REF1PHI ==> RF1PHI */
/* OM REF2PHI ==> RF2PHI */
/* OM NREF1PHI ==> NRF1PY */
/* OM NREF2PHI ==> NRF2PY */
/* OM DGDANGLE ==> DGDANG */
/* OM IRESTMAX ==> IRSTMX */
/* OM */
/* OM */
/* Subroutine */ int clschna_(x, y, z, atmind, bond, angle, newx, newy, newz, 
	omega, iter, entzus, pomega, maxdt, nangle, maxg)
real *x, *y, *z;
integer *atmind;
real *bond, *angle, *newx, *newy, *newz, *omega;
integer *iter;
logical *entzus;
real *pomega, *maxdt, *nangle, *maxg;
{
    /* Initialized data */

    static real rts[50] = { (float)0.,(float)0.,(float)0.,(float)0.,(float)0.,
	    (float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(
	    float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)
	    0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(
	    float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)
	    0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(
	    float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)
	    0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0. };
    static integer mxitrf = 50;
    static integer maxnw = 5;
    static integer mxmodt = 4;
    static integer mxitmn = 4;
    static integer bsmxit = 10;
    static real eps3 = (float).002;
    static integer itmax = 0;

    /* Format strings */
    static char fmt_100[] = "";
    static char fmt_300[] = "";
    static char fmt_200[] = "";
    static char fmt_500[] = "";
    static char fmt_600[] = "";
    static char fmt_700[] = "";
    static char fmt_1100[] = "";
    static char fmt_1200[] = "";
    static char fmt_1300[] = "";
    static char fmt_1430[] = "";
    static char fmt_1420[] = "";
    static char fmt_1605[] = "";
    static char fmt_1700[] = "";
    static char fmt_1670[] = "";
    static char fmt_1680[] = "";
    static char fmt_1660[] = "";
    static char fmt_9001[] = "(\0020Warning from CLSCHN -- Boundaries failed\
 to match.\002,\002 ILIM = \002,i1,\002 NP = \002,i1,\002 G(PHI1) =\002,4(1x\
,1pg14.7))";
    static char fmt_2000[] = "";
    static char fmt_9002[] = "(\002 RR,QR1,QR2 = \002,1p3g14.7)";
    static char fmt_9003[] = "(\002 LB,UB = \002,1p2g14.7)";
    static char fmt_2900[] = "";
    static char fmt_3300[] = "";
    static char fmt_2950[] = "";
    static char fmt_3000[] = "";
    static char fmt_9004[] = "(\002 PHI1: \002,4(1pg14.7,1x))";
    static char fmt_3400[] = "";
    static char fmt_9005[] = "(\002 Bound for eq. 28 for PHI1(\002,i1,\002\
,\002,i1,\002) is \002,1pg14.7)";
    static char fmt_9006[] = "(\002 Bound for eq. 34 for PHI1(\002,i1,\002\
,\002,i1,\002) is \002,1pg14.7)";
    static char fmt_9007[] = "(\002 APHI and BPHI for GET-PHI-RANGES are \
\002,1p2g14.7)";
    static char fmt_9008[] = "(\0020Error in CLSCHN -- Fell through \002,\
\002PHI CONDITIONAL\002)";
    static char fmt_9009[] = "(\002 Phi1 values out of order PHI1(1,\002,i1\
,\002) = \002,1pg14.7,\002 PHI1(2,\002,i1,\002) = \002,1pg14.7)";
    static char fmt_9010[] = "(\002 Eq. 28 bounds:\002,i2,\002:\002,1p4g14.7)"
	    ;
    static char fmt_9011[] = "(\002 Eq. 34 bounds:\002,i2,\002:\002,1p4g14.7)"
	    ;
    static char fmt_4500[] = "";
    static char fmt_9012[] = "(\0020KOUNT exceeded 10*BSMXIT.  KOUNT = \002,\
i12)";
    static char fmt_9013[] = "(1x,4(a,i3))";
    static char fmt_4700[] = "";
    static char fmt_4800[] = "";
    static char fmt_4900[] = "";
    static char fmt_5300[] = "";
    static char fmt_5700[] = "";
    static char fmt_5000[] = "";
    static char fmt_5100[] = "";
    static char fmt_5200[] = "";
    static char fmt_5400[] = "";
    static char fmt_5500[] = "";
    static char fmt_5600[] = "";
    static char fmt_5800[] = "";
    static char fmt_5900[] = "";
    static char fmt_6000[] = "";
    static char fmt_217[] = "(\0020Warning from CLSCHN -- An odd number of s\
olutions \002,\002were found.\002)";
    static char fmt_6500[] = "";
    static char fmt_6600[] = "";
    static char fmt_6700[] = "";
    static char fmt_6900[] = "";
    static char fmt_7100[] = "";
    static char fmt_7202[] = "";
    static char fmt_7204[] = "";
    static char fmt_9014[] = "(\002 Search for G, Angle sign = \002,9i3,\002\
 G =\002,1pg14.7)";
    static char fmt_7400[] = "";
    static char fmt_7600[] = "";
    static char fmt_7500[] = "";
    static char fmt_7720[] = "";
    static char fmt_7740[] = "";
    static char fmt_8300[] = "";
    static char fmt_8100[] = "";
    static char fmt_9015[] = "(\002 G found for NANGLE(\002,i1,\002,\002,i1\
,\002)=\002,f7.3)";
    static char fmt_8200[] = "";
    static char fmt_8420[] = "";
    static char fmt_8440[] = "";
    static char fmt_8480[] = "";
    static char fmt_8500[] = "";
    static char fmt_9016[] = "(\002 DG/DANGLE =\002,3(/1x,1p3g14.7))";
    static char fmt_8900[] = "";
    static char fmt_9017[] = "(\002 G(\002,f7.4,\002) =\002,1pg14.7)";
    static char fmt_9018[] = "(\002 Shift for ANGLE(\002,i1,\002,\002,i1,\
\002)=\002,1pg14.7)";
    static char fmt_9500[] = "";
    static char fmt_9600[] = "";
    static char fmt_9700[] = "";
    static char fmt_11200[] = "";
    static char fmt_9019[] = "(\002Muller's method is being used to solve th\
e eqn.\002)";
    static char fmt_10000[] = "";
    static char fmt_9020[] = "(\002 Regula falsi being used to solve the equ\
ation.\002)";
    static char fmt_11300[] = "";
    static char fmt_10100[] = "";
    static char fmt_10300[] = "";
    static char fmt_10400[] = "";
    static char fmt_10500[] = "";
    static char fmt_10600[] = "";
    static char fmt_10700[] = "";
    static char fmt_10800[] = "";
    static char fmt_11000[] = "";
    static char fmt_9021[] = "(\002 Incrementing MAXIT to \002,i4,\002 as mo\
re roots are likely.\002)";
    static char fmt_11100[] = "";
    static char fmt_9022[] = "(\002 A Root is \002,1pg14.7,\002 F is \002,1p\
g14.7,i5,\002 iterations\002)";
    static char fmt_11700[] = "";
    static char fmt_9023[] = "(\002 A solution in unconstrained space is \
\002,1pg14.7,\002  IER = \002,i3)";
    static char fmt_11820[] = "";
    static char fmt_12100[] = "";
    static char fmt_12200[] = "";
    static char fmt_12300[] = "";
    static char fmt_9024[] = "(\002 Sign change indicates G has zero value\
.\002)";
    static char fmt_12400[] = "";
    static char fmt_12500[] = "";
    static char fmt_12600[] = "";
    static char fmt_12700[] = "";
    static char fmt_12800[] = "";
    static char fmt_13000[] = "";
    static char fmt_13100[] = "";
    static char fmt_9025[] = "(\002 A Minimum G is at \002,1pg14.7,\002 G \
is \002,1pg14.7,i5,\002 iterations\002)";
    static char fmt_13500[] = "";
    static char fmt_13700[] = "";
    static char fmt_9026[] = "(1x,1p4g14.7)";
    static char fmt_14000[] = "";
    static char fmt_14200[] = "";
    static char fmt_14300[] = "";
    static char fmt_14400[] = "";
    static char fmt_14500[] = "";
    static char fmt_14600[] = "";
    static char fmt_14700[] = "";
    static char fmt_15700[] = "";
    static char fmt_15800[] = "";
    static char fmt_15900[] = "";
    static char fmt_16000[] = "";
    static char fmt_16400[] = "";
    static char fmt_9027[] = "(1x,1pg14.7,\002should == COS(THETA(5)) = \002\
,1pg14.7)";
    static char fmt_16600[] = "";
    static char fmt_9028[] = "(\002 The following number should be zero \002\
,1pg14.7)";

    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi(), s_wsli(), do_lio(), e_wsli();
    double cos(), sin(), sqrt(), atan2(), asin(), r_sign();
    integer i_dnnt();

    /* Local variables */
    static integer i99987, i99985, i99981, i99979, i99977, i99969, i99975, 
	    i99956, i99628, i99626, i99624, i99612, i99608;
    static real delf, phi1[4]	/* was [2][2] */, aphi, bphi, oldg, gphi;
    static logical done;
    static integer ilim;
    static real sfrf, tref;
    static integer sgns[2];
    static real swrf, diff1, diff2, frac1, frac2;
    static logical fnd0g;
    static real disc1, disc2, gphi1[36]	/* was [2][2][3][3] */;
    static integer nphi1;
    static real oldp1[4]	/* was [2][2] */, omega1neg;
    static integer savi1, savi2, savi3, savi4, savi5, savi6, savi7, savi8, 
	    savi9;
    static logical retry_clschn__;
    static real a, b, c, g, h;
    static integer i, j;
    static real omega1pos, diff12, diff13, diff14, p[18]	/* was [3][6] 
	    */, q[9]	/* was [3][3] */, r[3], s[3], t[3], u[3];
    static integer debug;
    static real v[3], delta, w, diff23, sigma[3], hiphi, denom, diff24, 
	    diff34, theta[6], lensd, tempa[3];
    static logical abort;
    static integer refnp;
    static real prodg;
    static integer npold;
    static real ratio, fxdfx;
    static logical found;
    static real lenmv;
    static integer maxit;
    static real frtdefneg;
    static logical swap34;
    static integer irest, npnew;
    static real g1, g2;
    static integer i1, i2, i3, i4, i5, kount, nroot, i6, i7, i8, i9;
    static real u0[3], v0[3], w0[3], t1, t2, u6[3], v6[3], frtdefpos, p1frac, 
	    t3, t4, rqsum, q3[3], g3, g4, refom1, oldom1;
    extern /* Subroutine */ int getuv_();
    static integer npolh1;
    static real smlom1, rtemp1, lb, lambda, omprv1, omprv2, dgdang[9]	/* 
	    was [3][3] */, ta[3], dt, dw, dx[3];
    extern doublereal abmerf_();
    static integer np;
    static real fx, tb[3], comega[6], ub, difphi;
    static integer vardbg;
    static real rt, oldang[9]	/* was [3][3] */, ctheta[6], ghiphi, frtdef, 
	    delfpr, rf1phi1, rf2phi1, smlgee, savang, glophi, somega[6];
    static integer sgnref[2];
    static real dfprlm, phirng, stheta[6], dispsg;
    static integer iterct;
    static real rr, xx, prvang[9]	/* was [3][3] */;
    static integer modtri, itrphi;
    static real nrf1py1, nrf2py1;
    static integer nplitl, sgnsml[2];
    static real lowphi, gpprev, stepsd, phirts[50];
    static integer nsdstp, is1, is2, it1, it2;
    static real sg1, frtprv, sg3;
    static integer irstmx;
    static real frtdefprev, nmshft[9]	/* was [3][3] */, lnshft, rr2, f1r, 
	    f2r, f3r, f4r, qr1, qr2, yy, zz;
    static char buffer[200];
    static logical overlp, possen, negsen, muller, donerf, donenw, toobig, 
	    gbndok, moddun, ok, adthet, stopsd, drpstp, caryon;
    static real den, arf, brf;
    static integer ind;
    static real frf, grf;
    static integer ier;
    static real dfx;
    extern doublereal dot_();
    static real qr12, rho[3], qr22, frt, wrf, num, dxx, sqr;
    extern /* Subroutine */ int cprint_(), die_(), abmerfi_();
    static real ang1, ang2;
    static integer i99610, i99602, i99630, i99811, i99821, i99803, i99913, 
	    i99851, i99815, i99862, i99818, i99909, i99856, i99929, i99848, 
	    i99805, i99780, i99765, i99703, i99705, i99699, i99678, i99668, 
	    i99666, i99662, i99660, i99632, i99998, i99996, i99994, i99989;

    /* Fortran I/O blocks */
    static icilist io___37 = { 0, buffer, 0, "(2PG14.7)", 200, 1 };
    static icilist io___39 = { 0, buffer, 0, "(2PG14.7)", 200, 1 };
    static icilist io___54 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___55 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___56 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___57 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___58 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___59 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___60 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___61 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___62 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___63 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___64 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___65 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___66 = { 0, buffer, 0, fmt_9001, 200, 1 };
    static icilist io___90 = { 0, buffer, 0, fmt_9002, 200, 1 };
    static icilist io___106 = { 0, buffer, 0, fmt_9003, 200, 1 };
    static icilist io___117 = { 0, buffer, 0, fmt_9004, 200, 1 };
    static icilist io___118 = { 0, buffer, 0, "(A)", 200, 1 };
    static icilist io___124 = { 0, buffer, 0, fmt_9005, 200, 1 };
    static icilist io___125 = { 0, buffer, 0, fmt_9006, 200, 1 };
    static icilist io___126 = { 0, buffer, 0, fmt_9007, 200, 1 };
    static icilist io___127 = { 0, buffer, 0, fmt_9008, 200, 1 };
    static icilist io___128 = { 0, buffer, 0, fmt_9009, 200, 1 };
    static icilist io___138 = { 0, buffer, 0, fmt_9010, 200, 1 };
    static icilist io___140 = { 0, buffer, 0, fmt_9011, 200, 1 };
    static icilist io___145 = { 0, buffer, 0, fmt_9012, 200, 1 };
    static icilist io___150 = { 0, buffer, 0, fmt_9013, 200, 1 };
    static icilist io___155 = { 0, buffer, 0, fmt_217, 200, 1 };
    static icilist io___160 = { 0, buffer, 0, "(A)", 200, 1 };
    static icilist io___175 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___193 = { 0, buffer, 0, fmt_9014, 200, 1 };
    static icilist io___194 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___199 = { 0, buffer, 0, fmt_9015, 200, 1 };
    static icilist io___200 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___206 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___207 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___209 = { 0, buffer, 0, fmt_9016, 200, 1 };
    static icilist io___218 = { 0, buffer, 0, fmt_9017, 200, 1 };
    static icilist io___233 = { 0, buffer, 0, fmt_9018, 200, 1 };
    static icilist io___234 = { 0, buffer, 0, "(A,1PG12.5,0PF8.5)", 200, 1 };
    static icilist io___235 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___245 = { 0, buffer, 0, fmt_9019, 200, 1 };
    static icilist io___249 = { 0, buffer, 0, fmt_9020, 200, 1 };
    static icilist io___271 = { 0, buffer, 0, fmt_9021, 200, 1 };
    static icilist io___276 = { 0, buffer, 0, fmt_9022, 200, 1 };
    static icilist io___278 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___287 = { 0, buffer, 0, fmt_9023, 200, 1 };
    static icilist io___288 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___290 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___293 = { 0, buffer, 0, fmt_9024, 200, 1 };
    static icilist io___295 = { 0, buffer, 0, fmt_9025, 200, 1 };
    static icilist io___297 = { 0, buffer, 0, 0, 200, 1 };
    static icilist io___298 = { 0, buffer, 0, fmt_9026, 200, 1 };
    static icilist io___299 = { 0, buffer, 0, fmt_9026, 200, 1 };
    static icilist io___309 = { 0, buffer, 0, "(A,1PG14.7)", 200, 1 };
    static icilist io___311 = { 0, buffer, 0, "(A,1PG14.7)", 200, 1 };
    static icilist io___319 = { 0, buffer, 0, fmt_9027, 200, 1 };
    static icilist io___320 = { 0, buffer, 0, fmt_9028, 200, 1 };


    /* Assigned format variables */
    char *i99998_fmt, *i99994_fmt, *i99996_fmt, *i99989_fmt, *i99987_fmt, *
	    i99985_fmt, *i99981_fmt, *i99979_fmt, *i99977_fmt, *i99969_fmt, *
	    i99956_fmt, *i99929_fmt, *i99913_fmt, *i99909_fmt, *i99975_fmt, *
	    i99862_fmt, *i99851_fmt, *i99848_fmt, *i99818_fmt, *i99815_fmt, *
	    i99803_fmt, *i99805_fmt, *i99811_fmt, *i99780_fmt, *i99765_fmt, *
	    i99703_fmt, *i99705_fmt, *i99699_fmt, *i99678_fmt, *i99668_fmt, *
	    i99666_fmt, *i99662_fmt, *i99660_fmt, *i99632_fmt, *i99630_fmt, *
	    i99628_fmt, *i99626_fmt, *i99624_fmt, *i99612_fmt, *i99610_fmt, *
	    i99608_fmt, *i99602_fmt;

/* OM */
/* OM      IMPLICIT CHARACTER*1000(A-H,J-Z) */
/* OM */




/* OM OMLUPD BNJ */
















/* #include "values.inc" */
/* #include "dbg.inc" */
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
    nangle -= 4;
    --pomega;
    --entzus;
    --omega;
    --newz;
    --newy;
    --newx;
    angle -= 4;
    bond -= 4;
    atmind -= 4;
    --z;
    --y;
    --x;

    /* Function Body */



    debug = dbg_1.dbg_clschn__ % 10;
    vardbg = dbg_1.dbg_clschn__ / 10 % 10;

    if (*iter == 0) {
	i99998 = 0;
	i99998_fmt = fmt_100;
	i99994 = 0;
	i99994_fmt = fmt_300;
	goto L2700;
    }
L100:
    i99996 = 0;
    i99996_fmt = fmt_200;
    goto L4200;
L200:
    return 0;
L300:
    for (i = 1; i <= 3; ++i) {
	if (! entzus[i]) {
	    pomega[i] = (float)3.141592653589794;
	} else {
	    pomega[i] = (float)0.;
	}
/* L400: */
    }
    i99989 = 0;
    i99989_fmt = fmt_500;
    goto L800;
L500:
    i99987 = 0;
    i99987_fmt = fmt_600;
    goto L1000;
L600:
    adthet = FALSE_;
    smlgee = (float)10.;
    i99985 = 0;
    i99985_fmt = fmt_700;
    goto L1600;
L700:
    sgns[0] = 1;
    sgns[1] = 1;
    nroot = 0;
    np = 1;
    iterct = 0;
    itrphi = 0;
    switch ((int)i99998) {
	case 0: goto L100;
    }
L800:
    for (i1 = 1; i1 <= 3; ++i1) {
	for (i2 = 1; i2 <= 3; ++i2) {
	    nangle[i1 + i2 * 3] = angle[i1 + i2 * 3];
/* L850: */
	}
/* L900: */
    }
    switch ((int)i99989) {
	case 0: goto L500;
	case 1: goto L5000;
	case 2: goto L5400;
	case 3: goto L5800;
    }
L1000:
    i99981 = 0;
    i99981_fmt = fmt_1100;
    goto L2100;
L1100:
    i99979 = 0;
    i99979_fmt = fmt_1200;
    goto L2300;
L1200:
    i99977 = 0;
    i99977_fmt = fmt_1300;
    goto L2800;
L1300:
    switch ((int)i99987) {
	case 0: goto L600;
	case 1: goto L5100;
	case 2: goto L5500;
	case 3: goto L5900;
	case 4: goto L7202;
	case 5: goto L7400;
	case 6: goto L8100;
	case 7: goto L8420;
	case 8: goto L8480;
	case 9: goto L9500;
	case 10: goto L9700;
    }
L1400:
    dispsg = smlgee;
    smlgee = (float)0.;
    i__1 = nphi1;
    for (np = 1; np <= i__1; ++np) {
	for (is2 = 1; is2 <= 2; ++is2) {
	    for (is1 = 1; is1 <= 2; ++is1) {
		sgns[0] = 3 - (is1 << 1);
		sgns[1] = 3 - (is2 << 1);
		omega[1] = phi1[(np << 1) - 2];
L1410:
		if (omega[1] > phi1[(np << 1) - 1]) {
		    omega[1] = phi1[(np << 1) - 1];
		    i99969 = 0;
		    i99969_fmt = fmt_1430;
		} else {
		    i99969 = 1;
		    i99969_fmt = fmt_1420;
		}
		goto L14100;
L1420:
		s_wsfi(&io___37);
		d__1 = omega[1] / .017453292519943299;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&gphi, (ftnlen)sizeof(real));
		e_wsfi();
		cprint_(buffer, 200L);
		omega[1] += (float).02;
		goto L1410;
L1430:
		s_wsfi(&io___39);
		d__1 = omega[1] / .017453292519943299;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&gphi, (ftnlen)sizeof(real));
		e_wsfi();
		cprint_(buffer, 200L);
/* L1440: */
	    }
/* L1450: */
	}
/* L1500: */
    }
    smlgee = dispsg;
    switch ((int)i99975) {
	case 0: goto L3400;
    }
L1600:
    i__1 = nphi1;
    for (np = 1; np <= i__1; ++np) {
	for (is1 = 1; is1 <= 2; ++is1) {
	    sgns[0] = 3 - (is1 << 1);
	    for (is2 = 1; is2 <= 2; ++is2) {
		sgns[1] = 3 - (is2 << 1);
		for (ilim = 1; ilim <= 2; ++ilim) {
		    omega[1] = phi1[ilim + (np << 1) - 3];
		    i99969 = 2;
		    i99969_fmt = fmt_1605;
		    goto L14100;
L1605:
		    gphi1[ilim + (np + (sgns[0] + sgns[1] * 3 << 1) << 1) + 
			    13] = gphi;
/* L1610: */
		}
/* L1620: */
	    }
/* L1650: */
	}
	if (phi1[(np << 1) - 2] != -3.141592653589794 || phi1[(np << 1) - 1] 
		!= 3.141592653589794) {
	    for (ilim = 1; ilim <= 2; ++ilim) {
		g1 = gphi1[ilim + (np + 8 << 1) + 13];
		g2 = gphi1[ilim + (np + 4 << 1) + 13];
		g3 = gphi1[ilim + (np - 4 << 1) + 13];
		g4 = gphi1[ilim + (np - 8 << 1) + 13];
		diff12 = g1 - g2;
		diff13 = g1 - g3;
		diff14 = g1 - g4;
		diff23 = g2 - g3;
		diff24 = g2 - g4;
		diff34 = g3 - g4;
		if (dabs(diff12) > eps3) {
		    if (dabs(diff13) <= eps3) {
			if (dabs(diff24) <= eps3) {
			    goto L1670;
			}
			i99956 = 0;
			i99956_fmt = fmt_1670;
			goto L1900;
		    } else if (dabs(diff14) > eps3) {
			i99956 = 1;
			i99956_fmt = fmt_1700;
			goto L1900;
		    } else {
			if (dabs(diff23) <= eps3) {
			    goto L1680;
			}
			i99956 = 2;
			i99956_fmt = fmt_1680;
			goto L1900;
		    }
		} else if (dabs(diff34) > eps3) {
		    i99956 = 3;
		    i99956_fmt = fmt_1660;
		    goto L1900;
		}
L1660:
		if (g1 * g2 < (float)0.) {
		    if (g1 >= (float)0.) {
			if (debug > 0) {
			    s_wsli(&io___54);
			    do_lio(&c__9, &c__1, "G1 and G2 are being adjust\
ed down.", 34L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np + 8 << 1) + 13] = (float)
				1.9999999999999999e-6;
			gphi1[ilim + (np + 4 << 1) + 13] = (float)
				4.9999999999999998e-7;
		    } else {
			if (debug > 0) {
			    s_wsli(&io___55);
			    do_lio(&c__9, &c__1, "G1 and G2 are being adjust\
ed up.", 32L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np + 8 << 1) + 13] = (float)
				4.9999999999999998e-7;
			gphi1[ilim + (np + 4 << 1) + 13] = (float)
				1.9999999999999999e-6;
		    }
		}
		if (g3 * g4 < (float)0.) {
		    if (g3 >= (float)0.) {
			if (debug > 0) {
			    s_wsli(&io___56);
			    do_lio(&c__9, &c__1, "G3 and G4 are being adjust\
ed down.", 34L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np - 4 << 1) + 13] = (float)
				1.9999999999999999e-6;
			gphi1[ilim + (np - 8 << 1) + 13] = (float)
				4.9999999999999998e-7;
		    } else {
			if (debug > 0) {
			    s_wsli(&io___57);
			    do_lio(&c__9, &c__1, "G3 and G4 are being adjust\
ed up.", 32L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np - 4 << 1) + 13] = (float)
				4.9999999999999998e-7;
			gphi1[ilim + (np - 8 << 1) + 13] = (float)
				1.9999999999999999e-6;
		    }
		}
		goto L1700;
L1670:
		if (g1 * g3 < (float)0.) {
		    if (g1 >= (float)0.) {
			if (debug > 0) {
			    s_wsli(&io___58);
			    do_lio(&c__9, &c__1, "G1 and G3 are being adjust\
ed down.", 34L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np + 8 << 1) + 13] = (float)
				1.9999999999999999e-6;
			gphi1[ilim + (np - 4 << 1) + 13] = (float)
				4.9999999999999998e-7;
		    } else {
			if (debug > 0) {
			    s_wsli(&io___59);
			    do_lio(&c__9, &c__1, "G1 and G3 are being adjust\
ed up.", 32L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np + 8 << 1) + 13] = (float)
				4.9999999999999998e-7;
			gphi1[ilim + (np - 4 << 1) + 13] = (float)
				1.9999999999999999e-6;
		    }
		}
		if (g2 * g4 < (float)0.) {
		    if (g2 >= (float)0.) {
			if (debug > 0) {
			    s_wsli(&io___60);
			    do_lio(&c__9, &c__1, "G2 and G4 are being adjust\
ed down.", 34L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np + 4 << 1) + 13] = (float)
				1.9999999999999999e-6;
			gphi1[ilim + (np - 8 << 1) + 13] = (float)
				4.9999999999999998e-7;
		    } else {
			if (debug > 0) {
			    s_wsli(&io___61);
			    do_lio(&c__9, &c__1, "G2 and G4 are being adjust\
ed up.", 32L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np + 4 << 1) + 13] = (float)
				4.9999999999999998e-7;
			gphi1[ilim + (np - 8 << 1) + 13] = (float)
				1.9999999999999999e-6;
		    }
		}
		goto L1700;
L1680:
		if (g1 * g4 < (float)0.) {
		    if (g1 >= (float)0.) {
			if (debug > 0) {
			    s_wsli(&io___62);
			    do_lio(&c__9, &c__1, "G1 and G4 are being adjust\
ed down.", 34L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np + 8 << 1) + 13] = (float)
				1.9999999999999999e-6;
			gphi1[ilim + (np - 8 << 1) + 13] = (float)
				4.9999999999999998e-7;
		    } else {
			if (debug > 0) {
			    s_wsli(&io___63);
			    do_lio(&c__9, &c__1, "G1 and G4 are being adjust\
ed up.", 32L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np + 8 << 1) + 13] = (float)
				4.9999999999999998e-7;
			gphi1[ilim + (np - 8 << 1) + 13] = (float)
				1.9999999999999999e-6;
		    }
		}
		if (g3 * g2 < (float)0.) {
		    if (g3 >= (float)0.) {
			if (debug > 0) {
			    s_wsli(&io___64);
			    do_lio(&c__9, &c__1, "G2 and G3 are being adjust\
ed down.", 34L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np - 4 << 1) + 13] = (float)
				1.9999999999999999e-6;
			gphi1[ilim + (np + 4 << 1) + 13] = (float)
				4.9999999999999998e-7;
		    } else {
			if (debug > 0) {
			    s_wsli(&io___65);
			    do_lio(&c__9, &c__1, "G2 and G3 are being adjust\
ed up.", 32L);
			    e_wsli();
			    cprint_(buffer, 200L);
			}
			gphi1[ilim + (np - 4 << 1) + 13] = (float)
				4.9999999999999998e-7;
			gphi1[ilim + (np + 4 << 1) + 13] = (float)
				1.9999999999999999e-6;
		    }
		}
L1700:
		;
	    }
	}
/* L1800: */
    }
    switch ((int)i99985) {
	case 0: goto L700;
	case 1: goto L5200;
	case 2: goto L5600;
	case 3: goto L6000;
	case 4: goto L6700;
	case 5: goto L6900;
	case 6: goto L7500;
	case 7: goto L8200;
    }
L1900:
    s_wsfi(&io___66);
    do_fio(&c__1, (char *)&ilim, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&np, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&g1, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&g2, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&g3, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&g4, (ftnlen)sizeof(real));
    e_wsfi();
    cprint_(buffer, 200L);
    i99929 = 0;
    i99929_fmt = fmt_2000;
    goto L3500;
L2000:
    switch ((int)i99956) {
	case 0: goto L1670;
	case 1: goto L1700;
	case 2: goto L1680;
	case 3: goto L1660;
    }
L2100:
    for (i = 1; i <= 3; ++i) {
	p[((i << 1) - 2) * 3] = bond[i * 3 + 1] - bond[i * 3 + 2] * cos(
		nangle[i * 3 + 1]);
	p[((i << 1) - 2) * 3 + 1] = bond[i * 3 + 2] * sin(nangle[i * 3 + 1]);
	p[((i << 1) - 2) * 3 + 2] = (float)0.;
	p[((i << 1) - 1) * 3] = bond[i * 3 + 3];
	p[((i << 1) - 1) * 3 + 1] = (float)0.;
	p[((i << 1) - 1) * 3 + 2] = (float)0.;
/* L2200: */
    }
    switch ((int)i99981) {
	case 0: goto L1100;
    }
L2300:
    for (i = 1; i <= 3; ++i) {
	theta[(i << 1) - 1] = 3.141592653589794 - nangle[i * 3 + 3];
	if (pomega[i] != (float)0.) {
	    theta[(i << 1) - 2] = nangle[i * 3 + 2] - nangle[i * 3 + 1];
	} else {
	    theta[(i << 1) - 2] = 6.283185307179588 - nangle[i * 3 + 1] - 
		    nangle[i * 3 + 2];
	}
/* L2400: */
    }
    for (i = 0; i <= 5; ++i) {
	ctheta[i] = cos(theta[i]);
	stheta[i] = sin(theta[i]);
/* L2500: */
    }
    for (i = 0; i <= 2; ++i) {
	q[i * 3] = p[((i << 1) + 1) * 3] * ctheta[i * 2] + p[i * 6];
	q[i * 3 + 1] = p[((i << 1) + 1) * 3] * stheta[i * 2] + p[i * 6 + 1];
	q[i * 3 + 2] = (float)0.;
	rho[i] = q[i * 3];
	sigma[i] = q[i * 3 + 1];
/* L2600: */
    }
    switch ((int)i99979) {
	case 0: goto L1200;
    }
L2700:
    getuv_(&x[1], &y[1], &z[1], &atmind[4], u0, v0);
    w0[0] = u0[1] * v0[2] - u0[2] * v0[1];
    w0[1] = u0[2] * v0[0] - u0[0] * v0[2];
    w0[2] = u0[0] * v0[1] - u0[1] * v0[0];
    getuv_(&x[1], &y[1], &z[1], &atmind[22], u6, v6);
    dx[0] = x[atmind[22]] - x[atmind[4]];
    dx[1] = y[atmind[22]] - y[atmind[4]];
    dx[2] = z[atmind[22]] - z[atmind[4]];
    s[0] = dot_(dx, u0, &c__3);
    s[1] = dot_(dx, v0, &c__3);
    s[2] = dot_(dx, w0, &c__3);
    u[0] = dot_(u6, u0, &c__3);
    u[1] = dot_(u6, v0, &c__3);
    u[2] = dot_(u6, w0, &c__3);
    v[0] = dot_(v6, u0, &c__3);
    v[1] = dot_(v6, v0, &c__3);
    v[2] = dot_(v6, w0, &c__3);
    switch ((int)i99994) {
	case 0: goto L300;
    }
L2800:
/* Computing 2nd power */
    r__1 = s[0] - q[0];
/* Computing 2nd power */
    r__2 = s[1] - q[1];
/* Computing 2nd power */
    r__3 = s[2] - q[2];
    rr2 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
    rr = sqrt(rr2);
/* Computing 2nd power */
    r__1 = q[3];
/* Computing 2nd power */
    r__2 = q[4];
/* Computing 2nd power */
    r__3 = q[5];
    qr12 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
    qr1 = sqrt(qr12);
/* Computing 2nd power */
    r__1 = q[6];
/* Computing 2nd power */
    r__2 = q[7];
/* Computing 2nd power */
    r__3 = q[8];
    qr22 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
    qr2 = sqrt(qr22);
    if (debug > 0 || vardbg > 0) {
	s_wsfi(&io___90);
	do_fio(&c__1, (char *)&rr, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&qr1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&qr2, (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 200L);
    }
/* Computing 2nd power */
    r__1 = qr1 + qr2;
    disc1 = r__1 * r__1 - rr2;
/* Computing 2nd power */
    r__1 = qr1 - qr2;
    disc2 = rr2 - r__1 * r__1;
    if (disc1 >= (float)0. && disc2 >= (float)0.) {
	rqsum = rr2 + qr12 - qr22;
	denom = qr12 * 2;
	t1 = rho[1] * rqsum;
	t2 = sigma[1] * sqrt(disc1) * sqrt(disc2);
	f1r = (t1 + t2) / denom;
	f2r = (t1 - t2) / denom;
	swap34 = FALSE_;
	t1 = sigma[2] * stheta[3];
	if (t1 < (float)0.) {
	    swap34 = ! swap34;
	}
	t2 = rho[2] * ctheta[3];
	t3 = -(doublereal)rho[1] * ctheta[2] + rqsum / 2 / sigma[1] * stheta[
		2] - sigma[1] * stheta[2];
	denom = rho[1] * stheta[2] / sigma[1] - ctheta[2];
	if (denom < (float)0.) {
	    swap34 = ! swap34;
	}
	f3r = (t3 + t1 - t2) / denom;
	f4r = (t3 - t1 - t2) / denom;
	if (swap34) {
	    t4 = f3r;
	    f3r = f4r;
	    f4r = t4;
	}
	ub = dmin(f1r,f3r);
	lb = dmax(f2r,f4r);
	if (debug > 0 || vardbg > 0) {
	    s_wsfi(&io___106);
	    do_fio(&c__1, (char *)&lb, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ub, (ftnlen)sizeof(real));
	    e_wsfi();
	    cprint_(buffer, 200L);
	}
	if (lb <= ub) {
	    for (i = 1; i <= 3; ++i) {
		tempa[i - 1] = s[i - 1] - q[i - 1];
/* L2820: */
	    }
	    t[0] = ctheta[0] * tempa[0] + stheta[0] * tempa[1];
	    t[1] = -(doublereal)stheta[0] * tempa[0] + ctheta[0] * tempa[1];
	    t[2] = tempa[2];
/* Computing 2nd power */
	    r__1 = t[1];
/* Computing 2nd power */
	    r__2 = t[2];
	    c = sqrt(r__1 * r__1 + r__2 * r__2);
	    aphi = (ub - t[0] * ctheta[1]) / (stheta[1] * c);
	    bphi = (lb - t[0] * ctheta[1]) / (stheta[1] * c);
	    delta = atan2(t[1], t[2]);
	    i99913 = 0;
	    i99913_fmt = fmt_2900;
	    goto L3700;
	} else {
	    nphi1 = 0;
	    goto L3200;
	}
    } else {
	nphi1 = 0;
	goto L3300;
    }
L2900:
    i__1 = nphi1;
    for (np = 1; np <= i__1; ++np) {
	if (phi1[(np << 1) - 2] == -3.141592653589794 && phi1[(np << 1) - 1] 
		== 3.141592653589794) {
	    goto L3100;
	}
	phi1[(np << 1) - 2] -= delta;
	phi1[(np << 1) - 1] -= delta;
	ilim = 1;
	i99909 = 0;
	i99909_fmt = fmt_2950;
	goto L3900;
L2950:
	if (rf2phi1 >= phi1[(np << 1) - 2]) {
/* Computing MIN */
	    r__1 = phi1[(np << 1) - 2];
	    phi1[(np << 1) - 2] = dmin(r__1,rf1phi1);
	} else {
	    phi1[(np << 1) - 2] = rf2phi1;
	}
	ilim = 2;
	i99909 = 1;
	i99909_fmt = fmt_3000;
	goto L3900;
L3000:
	if (rf1phi1 <= phi1[(np << 1) - 1]) {
/* Computing MAX */
	    r__1 = phi1[(np << 1) - 1];
	    phi1[(np << 1) - 1] = dmax(r__1,rf2phi1);
	} else {
	    phi1[(np << 1) - 1] = rf1phi1;
	}
L3100:
	;
    }
    if (debug > 0) {
	s_wsfi(&io___117);
	i__1 = nphi1;
	for (np = 1; np <= i__1; ++np) {
	    do_fio(&c__1, (char *)&phi1[(np << 1) - 2], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&phi1[(np << 1) - 1], (ftnlen)sizeof(real));
	}
	e_wsfi();
	cprint_(buffer, 200L);
    }
    if (nphi1 == 2) {
	if (phi1[1] > phi1[2]) {
	    if (debug > 0) {
		s_wsfi(&io___118);
		do_fio(&c__1, " PHI1 overlap corrected.", 24L);
		e_wsfi();
		cprint_(buffer, 200L);
	    }
	    phi1[1] = phi1[3];
	    nphi1 = 1;
	}
    }
L3200:
    if (debug > 0) {
	i99929 = 1;
	i99929_fmt = fmt_3300;
	goto L3500;
    }
L3300:
    if (debug > 0) {
	i99975 = 0;
	i99975_fmt = fmt_3400;
	goto L1400;
    }
L3400:
    switch ((int)i99977) {
	case 0: goto L1300;
    }
L3500:
    i__1 = nphi1;
    for (i = 1; i <= i__1; ++i) {
	for (j = 1; j <= 2; ++j) {
	    comega[0] = cos(phi1[j + (i << 1) - 3]);
	    somega[0] = sin(phi1[j + (i << 1) - 3]);
	    xx = t[0] * ctheta[1] + t[1] * stheta[1] * comega[0] + t[2] * 
		    stheta[1] * somega[0];
	    w = (rqsum - rho[1] * 2 * xx) / (sigma[1] * 2);
	    s_wsfi(&io___124);
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	    r__1 = rr2 - xx * xx - w * w;
	    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    e_wsfi();
	    cprint_(buffer, 200L);
	    s_wsfi(&io___125);
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	    r__1 = (rho[2] * ctheta[3] - (xx - rho[1]) * ctheta[2] - (w - 
		    sigma[1]) * stheta[2]) / (sigma[2] * stheta[3]);
	    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    e_wsfi();
	    cprint_(buffer, 200L);
/* L3550: */
	}
/* L3600: */
    }
    switch ((int)i99929) {
	case 0: goto L2000;
	case 1: goto L3300;
    }
L3700:
    if (debug > 0 || vardbg > 0) {
	s_wsfi(&io___126);
	do_fio(&c__1, (char *)&aphi, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&bphi, (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    if (bphi > aphi) {
	nphi1 = 0;
    } else if (bphi > (float)1.) {
	nphi1 = 0;
    } else if (aphi < (float)-1.) {
	nphi1 = 0;
    } else if (aphi == bphi) {
	nphi1 = 1;
	phi1[0] = asin(aphi);
	phi1[1] = phi1[0];
    } else if (aphi > (float)1. && bphi < (float)-1.) {
	nphi1 = 1;
	phi1[0] = (float)-3.141592653589794;
	phi1[1] = (float)3.141592653589794;
    } else if (aphi > (float)1. && bphi <= (float)1.) {
	nphi1 = 1;
	phi1[0] = asin(bphi);
	phi1[1] = 3.141592653589794 - phi1[0];
    } else if (bphi < (float)-1. && aphi >= (float)-1.) {
	nphi1 = 1;
	phi1[1] = asin(aphi);
	phi1[0] = -3.141592653589794 - phi1[1];
    } else if (aphi > (float)1. || bphi < (float)-1.) {
	s_wsfi(&io___127);
	e_wsfi();
	cprint_(buffer, 200L);
	die_();
    } else {
	nphi1 = 2;
	if (bphi < (float)0.) {
	    phi1[2] = asin(bphi);
	    phi1[3] = asin(aphi);
	    phi1[0] = -3.141592653589794 - phi1[3];
	    phi1[1] = -3.141592653589794 - phi1[2];
	} else {
	    phi1[0] = asin(bphi);
	    phi1[1] = asin(aphi);
	    phi1[2] = 3.141592653589794 - phi1[1];
	    phi1[3] = 3.141592653589794 - phi1[0];
	}
    }
    i__1 = nphi1;
    for (np = 1; np <= i__1; ++np) {
	if (phi1[(np << 1) - 2] > phi1[(np << 1) - 1]) {
	    s_wsfi(&io___128);
	    do_fio(&c__1, (char *)&np, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&phi1[(np << 1) - 2], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&np, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&phi1[(np << 1) - 1], (ftnlen)sizeof(real));
	    e_wsfi();
	    cprint_(buffer, 200L);
	}
/* L3800: */
    }
    switch ((int)i99913) {
	case 0: goto L2900;
    }
L3900:
    rf1phi1 = phi1[ilim + (np << 1) - 3];
    kount = 0;
L4000:
    comega[0] = cos(rf1phi1);
    somega[0] = sin(rf1phi1);
    xx = t[0] * ctheta[1] + t[1] * stheta[1] * comega[0] + t[2] * stheta[1] * 
	    somega[0];
    w = (rqsum - rho[1] * 2 * xx) / (sigma[1] * 2);
    dxx = stheta[1] * (t[2] * comega[0] - t[1] * somega[0]);
    dw = -(doublereal)rho[1] / sigma[1] * dxx;
    fx = rr2 - xx * xx - w * w + (float)1e-5;
    dfx = (xx * dxx + w * dw) * -2;
    if (dfx != (float)0.) {
	fxdfx = fx / dfx;
    } else {
	fxdfx = (float)0.;
    }
    nrf1py1 = rf1phi1 - fxdfx;
    ++kount;
    toobig = dabs(fxdfx) > (float).01 || dabs(fx) > (float).001;
    donenw = kount >= maxnw || nrf1py1 - rf1phi1 == (float)0. || dabs(fxdfx) 
	    <= (float)1e-6 || toobig;
    if (debug > 0) {
	s_wsfi(&io___138);
	do_fio(&c__1, (char *)&kount, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nrf1py1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rf1phi1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&fx, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dfx, (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    if (! toobig) {
	rf1phi1 = nrf1py1;
    }
    if (! donenw) {
	goto L4000;
    }
    rf2phi1 = phi1[ilim + (np << 1) - 3];
    kount = 0;
L4100:
    comega[0] = cos(rf2phi1);
    somega[0] = sin(rf2phi1);
    xx = t[0] * ctheta[1] + t[1] * stheta[1] * comega[0] + t[2] * stheta[1] * 
	    somega[0];
    w = (rqsum - rho[1] * 2 * xx) / (sigma[1] * 2);
    dxx = stheta[1] * (t[2] * comega[0] - t[1] * somega[0]);
    dw = -(doublereal)rho[1] / sigma[1] * dxx;
    fx = (rho[2] * ctheta[3] - (xx - rho[1]) * ctheta[2] - (w - sigma[1]) * 
	    stheta[2]) / (sigma[2] * stheta[3]);
    dfx = (-(doublereal)ctheta[2] * dxx - stheta[2] * dw) / (sigma[2] * 
	    stheta[3]);
    if (fx < (float)0.) {
	dfx = -(doublereal)dfx;
	fx = -(doublereal)fx;
    }
    fx = fx - (float)1. - (float)1e-5;
    if (dfx != (float)0.) {
	fxdfx = fx / dfx;
    } else {
	fxdfx = (float)0.;
    }
    nrf2py1 = rf2phi1 - fxdfx;
    ++kount;
    toobig = dabs(fxdfx) > (float).01 || dabs(fx) > (float).001;
    donenw = kount >= maxnw || nrf2py1 - rf2phi1 == (float)0. || dabs(fxdfx) 
	    <= (float)1e-6 || toobig;
    if (debug > 0) {
	s_wsfi(&io___140);
	do_fio(&c__1, (char *)&kount, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nrf2py1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rf2phi1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&fx, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dfx, (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    if (! toobig) {
	rf2phi1 = nrf2py1;
    }
    if (! donenw) {
	goto L4100;
    }
    if (rf1phi1 > rf2phi1) {
	tref = rf1phi1;
	rf1phi1 = rf2phi1;
	rf2phi1 = tref;
    }
    switch ((int)i99909) {
	case 0: goto L2950;
	case 1: goto L3000;
    }
L4200:
    done = FALSE_;
    found = FALSE_;
    goto L4400;
L4300:
    if (done) {
	switch ((int)i99996) {
	    case 0: goto L200;
	}
	goto L6100;
    }
L4400:
    if (np > nphi1) {
	goto L4600;
    }
    i99862 = 0;
    i99862_fmt = fmt_4500;
    goto L9900;
L4500:
    if (kount > bsmxit * 10) {
	s_wsfi(&io___145);
	do_fio(&c__1, (char *)&kount, (ftnlen)sizeof(integer));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    found = dabs(frtdef) <= (float)1e-6 && dabs(frt) <= (float)1e-6;
L4600:
    if (found) {
	if (debug > 0) {
	    itmax = max(itmax,kount);
	    irstmx = max(irest,irstmx);
	    s_wsfi(&io___150);
	    do_fio(&c__1, "KOUNT=", 6L);
	    do_fio(&c__1, (char *)&kount, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " ITMAX=", 7L);
	    do_fio(&c__1, (char *)&itmax, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " IREST =", 8L);
	    do_fio(&c__1, (char *)&irest, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " IRSTMX =", 9L);
	    do_fio(&c__1, (char *)&irstmx, (ftnlen)sizeof(integer));
	    e_wsfi();
	    cprint_(buffer, 200L);
	}
	done = TRUE_;
	++(*iter);
	++iterct;
	++itrphi;
	i99856 = 0;
	goto L15600;
    } else if (sgns[0] == 1) {
	sgns[0] = -1;
	nroot = 0;
    } else if (sgns[1] == 1) {
	sgns[1] = -1;
	sgns[0] = 1;
	nroot = 0;
    } else if (np < nphi1) {
	i99851 = 0;
	i99851_fmt = fmt_4700;
	goto L6100;
    } else if (*maxdt > (float)0. && ! adthet && itrphi == 0) {
	i99848 = 0;
	i99848_fmt = fmt_4800;
	goto L6200;
    } else if (pomega[1] == (float)0.) {
	i99851 = 1;
	i99851_fmt = fmt_4900;
	goto L6100;
    } else if (pomega[2] == (float)0.) {
	i99851 = 2;
	i99851_fmt = fmt_5300;
	goto L6100;
    } else if (pomega[3] != (float)0.) {
	done = TRUE_;
	*iter = 0;
    } else {
	i99851 = 3;
	i99851_fmt = fmt_5700;
	goto L6100;
    }
    goto L4300;
L4700:
    ++np;
    sgns[0] = 1;
    sgns[1] = 1;
    nroot = 0;
    goto L4300;
L4800:
    adthet = TRUE_;
    if (retry_clschn__) {
	np = 1;
	sgns[0] = 1;
	sgns[1] = 1;
	nroot = 0;
    }
    goto L4300;
L4900:
    pomega[1] = (float)3.141592653589794;
    adthet = FALSE_;
    i99989 = 1;
    i99989_fmt = fmt_5000;
    goto L800;
L5000:
    i99987 = 1;
    i99987_fmt = fmt_5100;
    goto L1000;
L5100:
    smlgee = (float)10.;
    i99985 = 1;
    i99985_fmt = fmt_5200;
    goto L1600;
L5200:
    itrphi = 0;
    np = 1;
    sgns[0] = 1;
    sgns[1] = 1;
    nroot = 0;
    goto L4300;
L5300:
    pomega[2] = (float)3.141592653589794;
    if (entzus[1]) {
	pomega[1] = (float)0.;
    }
    adthet = FALSE_;
    i99989 = 2;
    i99989_fmt = fmt_5400;
    goto L800;
L5400:
    i99987 = 2;
    i99987_fmt = fmt_5500;
    goto L1000;
L5500:
    smlgee = (float)10.;
    i99985 = 2;
    i99985_fmt = fmt_5600;
    goto L1600;
L5600:
    itrphi = 0;
    np = 1;
    sgns[0] = 1;
    sgns[1] = 1;
    nroot = 0;
    goto L4300;
L5700:
    pomega[3] = (float)3.141592653589794;
    if (entzus[2]) {
	pomega[2] = (float)0.;
    }
    if (entzus[1]) {
	pomega[1] = (float)0.;
    }
    adthet = FALSE_;
    i99989 = 3;
    i99989_fmt = fmt_5800;
    goto L800;
L5800:
    i99987 = 3;
    i99987_fmt = fmt_5900;
    goto L1000;
L5900:
    smlgee = (float)10.;
    i99985 = 3;
    i99985_fmt = fmt_6000;
    goto L1600;
L6000:
    itrphi = 0;
    np = 1;
    sgns[0] = 1;
    sgns[1] = 1;
    nroot = 0;
    goto L4300;
/* OM  OMLUPD BNJ ditch the variable format doobry for now... */
/* OM */
    if (! (iterct % 2 != 0)) {
	goto L6100;
    }
    s_wsfi(&io___155);
    e_wsfi();
/* OM      WRITE (BUFFER,217) ITERCT,((((GPHI1(ILIM,NP,IS1,IS2),ILIM=1,2),
 */
/* OM     2                       IS1=1,-1,-2),IS2=1,-1,-2),NP=1,NPHI1) */
/* OM  217 FORMAT('0Warning from CLSCHN -- An odd number of solutions ', 
*/
/* OM     2       'were found. ITERCT = ',I3/' GPHI1 values :' */
/* OM     3       <4*NPHI1>(/2(1X,1PG14.7))) */
    cprint_(buffer, 200L);
L6100:
    iterct = 0;
    switch ((int)i99851) {
	case 0: goto L4700;
	case 1: goto L4900;
	case 2: goto L5300;
	case 3: goto L5700;
    }
L6200:
    gbndok = TRUE_;
    if (dabs(smlgee) > *maxg) {
	goto L6800;
    }
    modtri = 0;
/* Computing MAX */
    r__1 = (float).001, r__2 = *maxdt / (float)10.;
    dt = dmax(r__1,r__2);
    abort = FALSE_;
    if (dabs(smlgee) <= (float)1e-6) {
	s_wsfi(&io___160);
	do_fio(&c__1, "Warning from CLSCHN - Solution was missed.", 42L);
	e_wsfi();
	cprint_(buffer, 200L);
    }
    goto L6400;
L6300:
    if (moddun) {
	goto L6800;
    }
/* OM */
/* OM OMUPD DW 19/11/91 CONVERT TO MORE READABLE IF () THEN CONSTRUCT */
/* OM */
/* OM 99826 MODDUN = ABS(SMLGEE) .LE. EPS2 .OR. */
/* OM     2              MODTRI .GE. MXMODT .OR. */
/* OM     3              ABORT */
/* OM */
L6400:
    if (dabs(smlgee) <= (float)1e-6 || modtri >= mxmodt || abort) {
	moddun = TRUE_;
    } else {
	moddun = FALSE_;
    }
/* OM */
    if (! moddun) {
	++modtri;
	if (nphi1 != 0) {
	    sg1 = smlgee;
	    npolh1 = nphi1;
	    i__1 = nphi1;
	    for (np = 1; np <= i__1; ++np) {
		oldp1[(np << 1) - 2] = phi1[(np << 1) - 2];
		oldp1[(np << 1) - 1] = phi1[(np << 1) - 1];
/* L6420: */
	    }
	    npold = nplitl;
	    oldom1 = smlom1;
	    i99818 = 0;
	    i99818_fmt = fmt_6500;
	    goto L8400;
	} else {
	    i99821 = 0;
	    goto L7000;
	}
    }
    goto L6300;
L6500:
    if (dabs(smlgee) <= (float)1e-6) {
	goto L6300;
    }
    i99815 = 0;
    i99815_fmt = fmt_6600;
    goto L9000;
L6600:
    if (dabs(smlgee) <= (float)1e-6) {
	goto L6300;
    }
    i99985 = 4;
    i99985_fmt = fmt_6700;
    goto L1600;
L6700:
    i99811 = 0;
    goto L11800;
L6800:
    retry_clschn__ = dabs(smlgee) <= (float)1e-6;
    if (retry_clschn__ && ! gbndok) {
	i99985 = 5;
	i99985_fmt = fmt_6900;
	goto L1600;
    }
L6900:
    switch ((int)i99848) {
	case 0: goto L4800;
    }
L7000:
    if (*maxg <= (float)1e3) {
	i99803 = 0;
	i99803_fmt = fmt_7100;
	goto L7700;
    } else {
	i99805 = 0;
	i99805_fmt = fmt_7100;
	goto L7200;
    }
L7100:
    switch ((int)i99821) {
	case 0: goto L6300;
    }
L7200:
    if (vardbg > 0) {
	s_wsli(&io___175);
	do_lio(&c__9, &c__1, "Searching for G thoroughly", 26L);
	e_wsli();
	cprint_(buffer, 200L);
    }
    debug = -debug;
    vardbg = -vardbg;
    sg3 = (float)10.;
    for (i1 = -1; i1 <= 1; i1 += 2) {
	nangle[4] = angle[4] + i1 * *maxdt;
	for (i2 = -1; i2 <= 1; i2 += 2) {
	    nangle[5] = angle[5] + i2 * *maxdt;
	    for (i3 = -1; i3 <= 1; i3 += 2) {
		nangle[6] = angle[6] + i3 * *maxdt;
		for (i4 = -1; i4 <= 1; i4 += 2) {
		    nangle[7] = angle[7] + i4 * *maxdt;
		    for (i5 = -1; i5 <= 1; i5 += 2) {
			nangle[8] = angle[8] + i5 * *maxdt;
			for (i6 = -1; i6 <= 1; i6 += 2) {
			    nangle[9] = angle[9] + i6 * *maxdt;
			    for (i7 = -1; i7 <= 1; i7 += 2) {
				nangle[10] = angle[10] + i7 * *maxdt;
				for (i8 = -1; i8 <= 1; i8 += 2) {
				    nangle[11] = angle[11] + i8 * *maxdt;
				    for (i9 = -1; i9 <= 1; i9 += 2) {
					nangle[12] = angle[12] + i9 * *maxdt;
					i99987 = 4;
					i99987_fmt = fmt_7202;
					goto L1000;
L7202:
					gbndok = FALSE_;
					sgns[0] = 1;
					sgns[1] = 1;
					i__1 = nphi1;
					for (np = 1; np <= i__1; ++np) {
					    omega[1] = phi1[(np << 1) - 2] + (
						    phi1[(np << 1) - 1] - 
						    phi1[(np << 1) - 2]) / (
						    float)2.;
					    i99969 = 3;
					    i99969_fmt = fmt_7204;
					    goto L14100;
L7204:
					    if (dabs(sg3) > dabs(gphi)) {
			  sg3 = gphi;
			  savi1 = i1;
			  savi2 = i2;
			  savi3 = i3;
			  savi4 = i4;
			  savi5 = i5;
			  savi6 = i6;
			  savi7 = i7;
			  savi8 = i8;
			  savi9 = i9;
					    }
/* L7206: */
					}
/* L7208: */
				    }
/* L7210: */
				}
/* L7212: */
			    }
/* L7214: */
			}
/* L7215: */
		    }
/* L7220: */
		}
/* L7240: */
	    }
/* L7250: */
	}
/* L7300: */
    }
    debug = -debug;
    vardbg = -vardbg;
    if (sg3 != (float)10.) {
	if (vardbg > 0) {
	    s_wsfi(&io___193);
	    do_fio(&c__1, (char *)&savi1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&savi2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&savi3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&savi4, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&savi5, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&savi6, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&savi7, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&savi8, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&savi9, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&sg3, (ftnlen)sizeof(real));
	    e_wsfi();
	    cprint_(buffer, 200L);
	}
	nangle[4] = angle[4] + savi1 * *maxdt;
	nangle[5] = angle[5] + savi2 * *maxdt;
	nangle[6] = angle[6] + savi3 * *maxdt;
	nangle[7] = angle[7] + savi4 * *maxdt;
	nangle[8] = angle[8] + savi5 * *maxdt;
	nangle[9] = angle[9] + savi6 * *maxdt;
	nangle[10] = angle[10] + savi7 * *maxdt;
	nangle[11] = angle[11] + savi8 * *maxdt;
	nangle[12] = angle[12] + savi9 * *maxdt;
	i99987 = 5;
	i99987_fmt = fmt_7400;
	goto L1000;
    } else {
	abort = TRUE_;
	goto L7600;
    }
L7400:
    i99985 = 6;
    i99985_fmt = fmt_7500;
    goto L1600;
L7500:
    gbndok = TRUE_;
    i99811 = 1;
    i99811_fmt = fmt_7600;
    goto L11800;
L7600:
    switch ((int)i99805) {
	case 0: goto L7100;
    }
L7700:
    if (vardbg > 0) {
	s_wsli(&io___194);
	do_lio(&c__9, &c__1, "Searching for G quickly", 23L);
	e_wsli();
	cprint_(buffer, 200L);
    }
    for (it1 = 1; it1 <= 3; ++it1) {
	for (it2 = 1; it2 <= 3; ++it2) {
	    savang = nangle[it1 + it2 * 3];
	    nangle[it1 + it2 * 3] = angle[it1 + it2 * 3] - *maxdt;
	    i99780 = 0;
	    i99780_fmt = fmt_7720;
	    goto L8000;
L7720:
	    if (nphi1 > 0) {
		goto L7900;
	    }
	    nangle[it1 + it2 * 3] = angle[it1 + it2 * 3] + *maxdt;
	    i99780 = 1;
	    i99780_fmt = fmt_7740;
	    goto L8000;
L7740:
	    if (nphi1 > 0) {
		goto L7900;
	    }
/* L7750: */
	}
/* L7800: */
    }
    abort = TRUE_;
L7900:
    switch ((int)i99803) {
	case 0: goto L7100;
    }
L8000:
    if (savang == nangle[it1 + it2 * 3]) {
	goto L8300;
    }
    i99987 = 6;
    i99987_fmt = fmt_8100;
    goto L1000;
L8100:
    gbndok = FALSE_;
    if (nphi1 <= 0) {
	goto L8300;
    }
    if (vardbg > 0) {
	s_wsfi(&io___199);
	do_fio(&c__1, (char *)&it1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&it2, (ftnlen)sizeof(integer));
	d__1 = nangle[it1 + it2 * 3] / .017453292519943299;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    i99985 = 7;
    i99985_fmt = fmt_8200;
    goto L1600;
L8200:
    gbndok = TRUE_;
    i99811 = 2;
    i99811_fmt = fmt_8300;
    goto L11800;
L8300:
    switch ((int)i99780) {
	case 0: goto L7720;
	case 1: goto L7740;
    }
L8400:
    if (vardbg > 0) {
	s_wsli(&io___200);
	do_lio(&c__9, &c__1, "Computing secants.", 18L);
	e_wsli();
	cprint_(buffer, 200L);
    }
    for (it1 = 1; it1 <= 3; ++it1) {
	for (it2 = 1; it2 <= 3; ++it2) {
	    ang1 = nangle[it1 + it2 * 3];
	    ang2 = ang1 + dt;
	    if (ang2 <= angle[it1 + it2 * 3] + *maxdt) {
		nangle[it1 + it2 * 3] = ang2;
		i99987 = 7;
		i99987_fmt = fmt_8420;
		goto L1000;
	    } else {
		ok = FALSE_;
		goto L8460;
	    }
L8420:
	    gbndok = FALSE_;
	    ok = nphi1 > 0;
	    if (! ok) {
		goto L8460;
	    }
	    i99765 = 0;
	    i99765_fmt = fmt_8440;
	    goto L8800;
L8440:
	    fnd0g = dabs(smlgee) < (float)1e-6 || smlgee * sg1 < (float)0. && 
		    nphi1 == npolh1 && nplitl == npold;
	    if (! fnd0g) {
		goto L8460;
	    }
	    smlgee = (float)0.;
	    if (vardbg > 0) {
		s_wsli(&io___206);
		do_lio(&c__9, &c__1, "Secant computation stopped as solution\
 was found.", 49L);
		e_wsli();
		cprint_(buffer, 200L);
	    }
	    goto L8700;
L8460:
	    if (ok) {
		goto L8520;
	    }
	    ang2 = ang1 - dt;
	    if (ang2 >= angle[it1 + it2 * 3] - *maxdt) {
		nangle[it1 + it2 * 3] = ang2;
		i99987 = 8;
		i99987_fmt = fmt_8480;
		goto L1000;
	    } else {
		ok = FALSE_;
		goto L8520;
	    }
L8480:
	    gbndok = FALSE_;
	    ok = nphi1 > 0;
	    if (! ok) {
		goto L8520;
	    }
	    i99765 = 1;
	    i99765_fmt = fmt_8500;
	    goto L8800;
L8500:
	    fnd0g = dabs(smlgee) < (float)1e-6 || smlgee * sg1 < (float)0. && 
		    nphi1 == npolh1 && nplitl == npold;
	    if (! fnd0g) {
		goto L8520;
	    }
	    smlgee = (float)0.;
	    if (vardbg > 0) {
		s_wsli(&io___207);
		do_lio(&c__9, &c__1, "Secant computation stopped as solution\
 was found.", 49L);
		e_wsli();
		cprint_(buffer, 200L);
	    }
	    goto L8700;
L8520:
	    nangle[it1 + it2 * 3] = ang1;
	    if (! ok) {
		dgdang[it1 + it2 * 3 - 4] = (float)0.;
	    } else {
		dgdang[it1 + it2 * 3 - 4] = (dabs(smlgee) - dabs(sg1)) / (
			ang2 - ang1);
	    }
/* L8550: */
	}
/* L8600: */
    }
    if (vardbg > 0) {
	s_wsfi(&io___209);
	for (it2 = 1; it2 <= 3; ++it2) {
	    for (it1 = 1; it1 <= 3; ++it1) {
		do_fio(&c__1, (char *)&dgdang[it1 + it2 * 3 - 4], (ftnlen)
			sizeof(real));
	    }
	}
	e_wsfi();
	cprint_(buffer, 200L);
    }
L8700:
    switch ((int)i99818) {
	case 0: goto L6500;
    }
L8800:
    sgns[0] = sgnsml[0];
    sgns[1] = sgnsml[1];
    if (npolh1 != 2 || nphi1 != 2) {
	denom = oldp1[(npolh1 << 1) - 1] - oldp1[0];
	if (denom != (float)0.) {
	    p1frac = (oldom1 - oldp1[0]) / denom;
	} else {
	    p1frac = (float)0.;
	}
	phirng = phi1[(nphi1 << 1) - 1] - phi1[0];
	if (nphi1 != 1) {
	    frac1 = (phi1[1] - phi1[0]) / phirng;
	    frac2 = (phi1[2] - phi1[0]) / phirng;
	    diff1 = p1frac - frac1;
	    diff2 = frac2 - p1frac;
	    if (diff1 * diff2 <= (float)0.) {
		omega[1] = phi1[0] + phirng * p1frac;
	    } else if (diff1 >= diff2) {
		omega[1] = phi1[2];
	    } else {
		omega[1] = phi1[1];
	    }
	    if (omega[1] > phi1[1]) {
		npnew = 2;
	    } else {
		npnew = 1;
	    }
	} else {
	    omega[1] = phi1[0] + phirng * p1frac;
	    npnew = 1;
	}
    } else {
	p1frac = (oldom1 - oldp1[(npold << 1) - 2]) / (oldp1[(npold << 1) - 1]
		 - oldp1[(npold << 1) - 2]);
	phirng = phi1[(npold << 1) - 1] - phi1[(npold << 1) - 2];
	omega[1] = phi1[(npold << 1) - 2] + p1frac * phirng;
	npnew = npold;
    }
    i99969 = 4;
    i99969_fmt = fmt_8900;
    goto L14100;
L8900:
    if (vardbg > 0) {
	s_wsfi(&io___218);
	do_fio(&c__1, (char *)&omega[1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&gphi, (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    smlgee = gphi;
    nplitl = npnew;
    switch ((int)i99765) {
	case 0: goto L8440;
	case 1: goto L8500;
	case 2: goto L9600;
    }
L9000:
    den = (float)0.;
    for (it1 = 1; it1 <= 3; ++it1) {
	for (it2 = 1; it2 <= 3; ++it2) {
/* Computing 2nd power */
	    r__1 = dgdang[it1 + it2 * 3 - 4];
	    den += r__1 * r__1;
	    oldang[it1 + it2 * 3 - 4] = nangle[it1 + it2 * 3];
/* L9050: */
	}
/* L9100: */
    }
    if (den >= (float)1e-6) {
	smlgee = sg1;
	npolh1 = nphi1;
	npold = nplitl;
	lensd = sqrt(den);
	stepsd = (float)1.;
	nsdstp = 0;
	goto L9300;
    } else {
	abort = TRUE_;
	switch ((int)i99815) {
	    case 0: goto L6600;
	}
	goto L9900;
    }
L9200:
    if (stopsd) {
	abort = lenmv < (float)1e-6;
	switch ((int)i99815) {
	    case 0: goto L6600;
	}
	goto L9900;
    }
L9300:
    ++nsdstp;
    oldg = smlgee;
    sgnref[0] = sgnsml[0];
    sgnref[1] = sgnsml[1];
    refom1 = smlom1;
    refnp = nplitl;
    lnshft = (float)0.;
    lenmv = (float)0.;
    for (it1 = 1; it1 <= 3; ++it1) {
	for (it2 = 1; it2 <= 3; ++it2) {
	    prvang[it1 + it2 * 3 - 4] = nangle[it1 + it2 * 3];
	    nangle[it1 + it2 * 3] = oldang[it1 + it2 * 3 - 4] - stepsd * 
		    dgdang[it1 + it2 * 3 - 4] * dabs(sg1) / den;
	    if (nangle[it1 + it2 * 3] < angle[it1 + it2 * 3] - *maxdt) {
		nangle[it1 + it2 * 3] = angle[it1 + it2 * 3] - *maxdt;
	    } else if (nangle[it1 + it2 * 3] > angle[it1 + it2 * 3] + *maxdt) 
		    {
		nangle[it1 + it2 * 3] = angle[it1 + it2 * 3] + *maxdt;
	    }
	    nmshft[it1 + it2 * 3 - 4] = nangle[it1 + it2 * 3] - prvang[it1 + 
		    it2 * 3 - 4];
/* Computing 2nd power */
	    r__1 = nmshft[it1 + it2 * 3 - 4];
	    lnshft += r__1 * r__1;
/* Computing 2nd power */
	    r__1 = nangle[it1 + it2 * 3] - oldang[it1 + it2 * 3 - 4];
	    lenmv += r__1 * r__1;
	    if (vardbg > 0) {
		s_wsfi(&io___233);
		do_fio(&c__1, (char *)&it1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nmshft[it1 + it2 * 3 - 4], (ftnlen)
			sizeof(real));
		e_wsfi();
		cprint_(buffer, 200L);
	    }
/* L9350: */
	}
/* L9400: */
    }
    i99987 = 9;
    i99987_fmt = fmt_9500;
    goto L1000;
L9500:
    gbndok = FALSE_;
    if (nphi1 != 0) {
	i99765 = 2;
	i99765_fmt = fmt_9600;
	goto L8800;
    } else {
	smlgee = (float)10.;
    }
L9600:
    fnd0g = dabs(smlgee) < (float)1e-6 || smlgee * sg1 < (float)0. && nphi1 ==
	     npolh1 && nplitl == npold;
    if (! fnd0g) {
	lenmv = sqrt(lenmv);
	lnshft = sqrt(lnshft);
	stepsd *= (float)1.5;
	if (vardbg > 0) {
	    s_wsfi(&io___234);
	    do_fio(&c__1, " LNSHFT,LENMV =", 15L);
	    do_fio(&c__1, (char *)&lnshft, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&lenmv, (ftnlen)sizeof(real));
	    e_wsfi();
	    cprint_(buffer, 200L);
	}
	if (dabs(smlgee) <= dabs(oldg) || nsdstp <= 1) {
	    goto L9800;
	}
	for (it1 = 1; it1 <= 3; ++it1) {
	    for (it2 = 1; it2 <= 3; ++it2) {
		nangle[it1 + it2 * 3] = prvang[it1 + it2 * 3 - 4];
/* L9620: */
	    }
/* L9650: */
	}
	i99987 = 10;
	i99987_fmt = fmt_9700;
	goto L1000;
    } else {
	smlgee = (float)0.;
	if (vardbg > 0) {
	    s_wsli(&io___235);
	    do_lio(&c__9, &c__1, "SD stopping as solution found.", 30L);
	    e_wsli();
	    cprint_(buffer, 200L);
	}
	goto L9800;
    }
L9700:
    gbndok = FALSE_;
    smlgee = oldg;
    sgnsml[0] = sgnref[0];
    sgnsml[1] = sgnref[1];
    smlom1 = refom1;
    nplitl = refnp;
L9800:
    drpstp = (dabs(smlgee) >= dabs(oldg) || nphi1 == 0) && nsdstp == 1 && 
	    lnshft >= (float)1e-6;
    if (! drpstp) {
	stopsd = dabs(smlgee) < (float)1e-6 || dabs(smlgee) >= dabs(oldg) || 
		lnshft < (float)1e-6;
    } else {
	stepsd = (float).1;
	stopsd = FALSE_;
    }
    goto L9200;
L9900:
    muller = TRUE_;
    kount = 0;
    irest = 0;
    lowphi = phi1[(np << 1) - 2];
    hiphi = phi1[(np << 1) - 1];
    difphi = hiphi - lowphi;
    glophi = gphi1[(np + (sgns[0] + sgns[1] * 3 << 1) << 1) + 14];
    ghiphi = gphi1[(np + (sgns[0] + sgns[1] * 3 << 1) << 1) + 15];
    if (nroot == 0) {
	if (dabs(glophi) <= (float)1e-6) {
	    rt = (float)-10.;
	    nroot = 1;
	    omega[1] = lowphi;
	    frtdef = glophi;
	    frt = glophi;
	    if (TRUE_) {
		goto L11200;
	    }
	} else if (dabs(ghiphi) <= (float)1e-6) {
	    rt = (float)10.;
	    nroot = 1;
	    omega[1] = hiphi;
	    frtdef = ghiphi;
	    frt = ghiphi;
	    if (TRUE_) {
		goto L11200;
	    }
	}
    }
    ++nroot;
    prodg = glophi * ghiphi;
    if (prodg >= (float)0. || nroot != 1) {
	if (debug > 0) {
	    s_wsfi(&io___245);
	    e_wsfi();
	    cprint_(buffer, 200L);
	}
	possen = FALSE_;
	negsen = FALSE_;
	rt = (float)-10.;
	i99703 = 0;
	i99703_fmt = fmt_10000;
	goto L13400;
    } else {
	if (debug > 0) {
	    s_wsfi(&io___249);
	    e_wsfi();
	    cprint_(buffer, 200L);
	}
	frf = glophi;
	grf = ghiphi;
	arf = lowphi;
	brf = hiphi;
	i99705 = 0;
	i99705_fmt = fmt_11300;
	goto L11500;
    }
    goto L11300;
L10000:
    rt = (float)10.;
    i99703 = 1;
    i99703_fmt = fmt_10100;
    goto L13400;
L10100:
    maxit = bsmxit;
    rts[nroot - 1] = (float)0.;
L10200:
    h = (float).5;
    rt = rts[nroot - 1] + h;
    i99703 = 2;
    i99703_fmt = fmt_10300;
    goto L13400;
L10300:
    i99699 = 0;
    i99699_fmt = fmt_10400;
    goto L11400;
L10400:
    delfpr = frtdef;
    rt = rts[nroot - 1] - h;
    i99703 = 3;
    i99703_fmt = fmt_10500;
    goto L13400;
L10500:
    i99699 = 1;
    i99699_fmt = fmt_10600;
    goto L11400;
L10600:
    frtprv = frtdef;
    delfpr = frtprv - delfpr;
    rt = rts[nroot - 1];
    i99703 = 4;
    i99703_fmt = fmt_10700;
    goto L13400;
L10700:
    i99699 = 2;
    i99699_fmt = fmt_10800;
    goto L11400;
L10800:
    lambda = (float)-.5;
L10900:
    delf = frtdef - frtprv;
    dfprlm = delfpr * lambda;
    num = -(doublereal)frtdef * (lambda + (float)1.) * (float)2.;
    g = (lambda * (float)2. + (float)1.) * delf - lambda * dfprlm;
    sqr = g * g + num * (float)2. * lambda * (delf - dfprlm);
    if (sqr < (float)0.) {
	sqr = (float)0.;
    }
    sqr = sqrt(sqr);
    den = g + sqr;
    if (g * sqr < (float)0.) {
	den = g - sqr;
    }
    if (dabs(den) == (float)0.) {
	den = (float)1.;
    }
    lambda = num / den;
    frtprv = frtdef;
    delfpr = delf;
    if (dabs(lambda) > (float)100.) {
	lambda = r_sign(&c_b352, &lambda);
    }
    h *= lambda;
    rt += h;
    if (kount <= maxit) {
	i99703 = 5;
	i99703_fmt = fmt_11000;
	goto L13400;
    } else if (kount <= bsmxit * 10) {
	caryon = possen && negsen || nroot % 2 == 0 && prodg > (float)0. || 
		nroot % 2 == 1 && prodg < (float)0.;
	if (possen && negsen) {
	    arf = omega1pos;
	    brf = omega1neg;
	    frf = frtdefpos;
	    grf = frtdefneg;
	    i99705 = 1;
	    i99705_fmt = fmt_11200;
	    goto L11500;
	} else if (caryon) {
	    if (debug > 0) {
		s_wsfi(&io___271);
		do_fio(&c__1, (char *)&maxit, (ftnlen)sizeof(integer));
		e_wsfi();
		cprint_(buffer, 200L);
	    }
	    maxit += bsmxit;
	    i99703 = 5;
	    i99703_fmt = fmt_11000;
	    goto L13400;
	}
    }
    goto L11200;
L11000:
/* Computing MAX */
    r__1 = dabs(frt), r__2 = dabs(frtdef);
    if (dmax(r__1,r__2) <= (float)1e-6) {
	goto L11200;
    }
    i99699 = 3;
    i99699_fmt = fmt_11100;
    goto L11400;
L11100:
/* Computing MAX */
    r__3 = (r__1 = omega[1] - omprv1, dabs(r__1)), r__4 = (r__2 = omprv1 - 
	    omprv2, dabs(r__2));
    overlp = dmax(r__3,r__4) == (float)0.;
    if (! overlp) {
	if (dabs(frtdef) <= dabs(frtprv) * (float)10.) {
	    goto L10900;
	}
	h /= (float)2.;
	lambda /= (float)2.;
	rt -= h;
	if (TRUE_) {
	    i99703 = 5;
	    goto L13400;
	}
    } else {
	if (rts[nroot - 1] <= (float)0.) {
	    rts[nroot - 1] = (float).6 - rts[nroot - 1];
	} else {
	    rts[nroot - 1] = -(doublereal)rts[nroot - 1];
	}
	++irest;
	goto L10200;
    }
L11200:
    rts[nroot - 1] = rt;
    phirts[nroot - 1] = omega[1];
    if (debug > 0) {
	s_wsfi(&io___276);
	do_fio(&c__1, (char *)&omega[1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&frt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&kount, (ftnlen)sizeof(integer));
	e_wsfi();
	cprint_(buffer, 200L);
    }
L11300:
    switch ((int)i99862) {
	case 0: goto L4500;
    }
L11400:
    if (frtdefprev * frtdef >= (float)0. || kount <= bsmxit) {
	switch ((int)i99699) {
	    case 0: goto L10400;
	    case 1: goto L10600;
	    case 2: goto L10800;
	    case 3: goto L11100;
	}
    } else {
	if (debug > 0) {
	    s_wsli(&io___278);
	    do_lio(&c__9, &c__1, "Switching to modified regula falsi.", 35L);
	    e_wsli();
	    cprint_(buffer, 200L);
	}
	frf = frtdef;
	grf = frtdefprev;
	arf = omega[1];
	brf = omprv1;
	i99705 = 1;
    }
L11500:
    muller = FALSE_;
    maxit = kount + mxitrf;
    wrf = arf;
    frtdef = frf;
    sfrf = r_sign(&c_b366, &frf);
L11600:
    wrf = (frf * brf - grf * arf) / (frf - grf);
    swrf = r_sign(&c_b366, &frtdef);
    omega[1] = wrf;
    rt = wrf;
    i99678 = 0;
    i99678_fmt = fmt_11700;
    goto L13600;
L11700:
    if (sfrf * frtdef >= (float)0.) {
	arf = wrf;
	frf = frtdef;
	if (swrf * frtdef > (float)0.) {
	    grf /= (float)2.;
	}
    } else {
	brf = wrf;
	grf = frtdef;
	if (swrf * frtdef > (float)0.) {
	    frf /= (float)2.;
	}
    }
/* Computing MAX */
    r__2 = dabs(frt), r__3 = dabs(frtdef);
    donerf = kount > maxit || dmax(r__2,r__3) <= (float)1e-6 || (r__1 = brf - 
	    arf, dabs(r__1)) <= (float)1e-6;
    if (! donerf) {
	goto L11600;
    }
    phirts[nroot - 1] = rt;
    if (dabs(frt) >= (float)1e-6) {
	ratio = frtdef / frt;
    } else {
	ratio = frtdef / (float)1e-6;
    }
    if ((r__1 = brf - arf, dabs(r__1)) <= (float)1e-6 && dabs(ratio) < (float)
	    100.) {
	frtdef = (float)1e-6;
	frt = (float)1e-6;
    }


/* OM OMUPD BNJ update code for erfi function... */

    rtemp1 = (rt - lowphi) * (float)2. / (hiphi - lowphi) - (float)1.;
    abmerfi_(&rtemp1, &rts[nroot - 1], &ier);

    if (debug > 1) {
	s_wsfi(&io___287);
	do_fio(&c__1, (char *)&rts[nroot - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    switch ((int)i99705) {
	case 0: goto L11300;
	case 1: goto L11200;
    }
L11800:
    if (vardbg > 0) {
	s_wsli(&io___288);
	do_lio(&c__9, &c__1, "Finding minimum G", 17L);
	e_wsli();
	cprint_(buffer, 200L);
    }
    smlgee = (float)10.;
    i__1 = nphi1;
    for (np = 1; np <= i__1; ++np) {
	for (is2 = 1; is2 >= -1; is2 += -2) {
	    sgns[1] = is2;
	    for (is1 = 1; is1 >= -1; is1 += -2) {
		sgns[0] = is1;
		i99668 = 0;
		i99668_fmt = fmt_11820;
		goto L12000;
L11820:
		;
	    }
/* L11850: */
	}
/* L11900: */
    }
    if (vardbg > 0) {
	s_wsli(&io___290);
	do_lio(&c__9, &c__1, "Smallest G = ", 13L);
	do_lio(&c__4, &c__1, (char *)&smlgee, (ftnlen)sizeof(real));
	e_wsli();
	cprint_(buffer, 200L);
    }
    switch ((int)i99811) {
	case 0: goto L6300;
	case 1: goto L7600;
	case 2: goto L8300;
    }
L12000:
    muller = TRUE_;
    kount = 0;
    lowphi = phi1[(np << 1) - 2];
    hiphi = phi1[(np << 1) - 1];
    difphi = hiphi - lowphi;
    glophi = gphi1[(np + (sgns[0] + sgns[1] * 3 << 1) << 1) + 14];
    ghiphi = gphi1[(np + (sgns[0] + sgns[1] * 3 << 1) << 1) + 15];
    gphi = glophi;
    omega[1] = lowphi;
    i99666 = 0;
    i99666_fmt = fmt_12100;
    goto L14800;
L12100:
    gphi = ghiphi;
    omega[1] = hiphi;
    i99666 = 1;
    i99666_fmt = fmt_12200;
    goto L14800;
L12200:
    if (dabs(smlgee) < (float)1e-6) {
	goto L13200;
    }
    nroot = 1;
    prodg = glophi * ghiphi;
    if (prodg >= (float)0.) {
	possen = FALSE_;
	negsen = FALSE_;
	maxit = mxitmn;
	rts[nroot - 1] = (float)0.;
	h = (float).5;
	rt = rts[nroot - 1] + h;
	i99662 = 0;
	i99662_fmt = fmt_12300;
	goto L13900;
    } else {
	if (vardbg > 0) {
	    s_wsfi(&io___293);
	    e_wsfi();
	    cprint_(buffer, 200L);
	}
	smlgee = (float)0.;
	goto L13200;
    }
L12300:
    i99660 = 0;
    i99660_fmt = fmt_12400;
    goto L13300;
L12400:
    delfpr = gphi;
    rt = rts[nroot - 1] - h;
    i99662 = 1;
    i99662_fmt = fmt_12500;
    goto L13900;
L12500:
    i99660 = 1;
    i99660_fmt = fmt_12600;
    goto L13300;
L12600:
    frtprv = gphi;
    delfpr = frtprv - delfpr;
    rt = rts[nroot - 1];
    i99662 = 2;
    i99662_fmt = fmt_12700;
    goto L13900;
L12700:
    i99660 = 2;
    i99660_fmt = fmt_12800;
    goto L13300;
L12800:
    lambda = (float)-.5;
L12900:
    delf = gphi - frtprv;
    dfprlm = delfpr * lambda;
    num = -(doublereal)gphi * (lambda + (float)1.) * (float)2.;
    g = (lambda * (float)2. + (float)1.) * delf - lambda * dfprlm;
    sqr = g * g + num * (float)2. * lambda * (delf - dfprlm);
    if (sqr < (float)0.) {
	sqr = (float)0.;
    }
    sqr = sqrt(sqr);
    den = g + sqr;
    if (g * sqr < (float)0.) {
	den = g - sqr;
    }
    if (dabs(den) == (float)0.) {
	den = (float)1.;
    }
    lambda = num / den;
    frtprv = gphi;
    delfpr = delf;
    if (dabs(lambda) > (float)100.) {
	lambda = r_sign(&c_b352, &lambda);
    }
    h *= lambda;
    rt += h;
    if (kount > maxit) {
	goto L13200;
    }
    i99662 = 3;
    i99662_fmt = fmt_13000;
    goto L13900;
L13000:
    if (dabs(smlgee) <= (float)1e-6) {
	goto L13200;
    }
    i99660 = 3;
    i99660_fmt = fmt_13100;
    goto L13300;
L13100:
/* Computing MAX */
    r__3 = (r__1 = omega[1] - omprv1, dabs(r__1)), r__4 = (r__2 = omprv1 - 
	    omprv2, dabs(r__2));
    overlp = dmax(r__3,r__4) == (float)0.;
    if (! overlp) {
	if (dabs(gphi) <= dabs(frtprv) * (float)10.) {
	    goto L12900;
	}
	h /= (float)2.;
	lambda /= (float)2.;
	rt -= h;
	if (TRUE_) {
	    i99662 = 3;
	    goto L13900;
	}
    }
L13200:
    if (vardbg > 0) {
	s_wsfi(&io___295);
	do_fio(&c__1, (char *)&omega[1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&smlgee, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&kount, (ftnlen)sizeof(integer));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    switch ((int)i99668) {
	case 0: goto L11820;
    }
L13300:
    if (gpprev * gphi >= (float)0.) {
	switch ((int)i99660) {
	    case 0: goto L12400;
	    case 1: goto L12600;
	    case 2: goto L12800;
	    case 3: goto L13100;
	}
    } else {
	if (vardbg > 0) {
	    s_wsli(&io___297);
	    do_lio(&c__9, &c__1, "Sign change seen so zero inferred.", 34L);
	    e_wsli();
	    cprint_(buffer, 200L);
	}
	smlgee = (float)0.;
	goto L13200;
    }
L13400:
    omprv2 = omprv1;
    omprv1 = omega[1];
    omega[1] = lowphi + difphi * (float).5 * (abmerf_(&rt) + (float)1.);
    frtdefprev = frtdef;
    i99678 = 1;
    i99678_fmt = fmt_13500;
    goto L13600;
L13500:
    if (frtdef > (float)0.) {
	possen = TRUE_;
	frtdefpos = frtdef;
	omega1pos = omega[1];
    } else if (frtdef < (float)0.) {
	negsen = TRUE_;
	frtdefneg = frtdef;
	omega1neg = omega[1];
    }
    switch ((int)i99703) {
	case 0: goto L10000;
	case 1: goto L10100;
	case 2: goto L10300;
	case 3: goto L10500;
	case 4: goto L10700;
	case 5: goto L11000;
    }
L13600:
    ++kount;
    if (omega[1] == lowphi) {
	gphi = glophi;
    } else if (omega[1] != hiphi) {
	i99969 = 5;
	i99969_fmt = fmt_13700;
	goto L14100;
    } else {
	gphi = ghiphi;
    }
L13700:
    frt = gphi;
    frtdef = frt;
    i__1 = nroot - 1;
    for (j = 1; j <= i__1; ++j) {
	den = omega[1] - phirts[j - 1];
	if (dabs(den) > (float)1e-6) {
	    frtdef /= den;
	} else if (! muller) {
	    den = (float)1e-6;
	    frtdef /= den;
	} else {
	    if (dabs(rt) <= (float)3.2) {
		rts[nroot - 1] = rt + (float).001;
	    } else {
		rts[nroot - 1] = ((integer) rt - rt) * 2;
	    }
	    if (debug > 1) {
		s_wsfi(&io___298);
		do_fio(&c__1, (char *)&rt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&frt, (ftnlen)sizeof(real));
		e_wsfi();
		cprint_(buffer, 200L);
	    }
	    if ((real) kount <= bsmxit * (float)11.) {
		goto L10200;
	    }
	    goto L11200;
	}
/* L13800: */
    }
    if (debug > 1 || vardbg > 1) {
	s_wsfi(&io___299);
	do_fio(&c__1, (char *)&rt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&omega[1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&frt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&frtdef, (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    switch ((int)i99678) {
	case 0: goto L11700;
	case 1: goto L13500;
    }
L13900:
    ++kount;
    gpprev = gphi;
    omega[1] = lowphi + difphi * (float).5 * (abmerf_(&rt) + (float)1.);
    if (omega[1] == lowphi) {
	gphi = glophi;
    } else if (omega[1] != hiphi) {
	i99969 = 6;
	i99969_fmt = fmt_14000;
	goto L14100;
    } else {
	gphi = ghiphi;
    }
L14000:
    switch ((int)i99662) {
	case 0: goto L12300;
	case 1: goto L12500;
	case 2: goto L12700;
	case 3: goto L13000;
    }
L14100:
    i99632 = 0;
    i99632_fmt = fmt_14200;
    goto L14900;
L14200:
    i99630 = 0;
    i99630_fmt = fmt_14300;
    goto L15100;
L14300:
    i99628 = 0;
    i99628_fmt = fmt_14400;
    goto L15200;
L14400:
    i99626 = 0;
    i99626_fmt = fmt_14500;
    goto L15300;
L14500:
    i99624 = 0;
    i99624_fmt = fmt_14600;
    goto L15400;
L14600:
    i99666 = 2;
    i99666_fmt = fmt_14700;
    goto L14800;
L14700:
    switch ((int)i99969) {
	case 0: goto L1430;
	case 1: goto L1420;
	case 2: goto L1605;
	case 3: goto L7204;
	case 4: goto L8900;
	case 5: goto L13700;
	case 6: goto L14000;
	case 7: goto L15700;
    }
L14800:
    if (dabs(gphi) < dabs(smlgee)) {
	smlgee = gphi;
	smlom1 = omega[1];
	sgnsml[0] = sgns[0];
	sgnsml[1] = sgns[1];
	nplitl = np;
    }
    switch ((int)i99666) {
	case 0: goto L12100;
	case 1: goto L12200;
	case 2: goto L14700;
    }
L14900:
    for (i = 1; i <= 3; ++i) {
	ta[i - 1] = s[i - 1] - q[i - 1];
/* L15000: */
    }
    r[0] = ctheta[0] * ta[0] + stheta[0] * ta[1];
    r[1] = -(doublereal)stheta[0] * ta[0] + ctheta[0] * ta[1];
    r[2] = ta[2];
    comega[0] = cos(omega[1]);
    somega[0] = sin(omega[1]);
    ta[0] = r[0];
    ta[1] = comega[0] * r[1] + somega[0] * r[2];
    ta[2] = -(doublereal)somega[0] * r[1] + comega[0] * r[2];
    r[0] = ctheta[1] * ta[0] + stheta[1] * ta[1];
    r[1] = -(doublereal)stheta[1] * ta[0] + ctheta[1] * ta[1];
    r[2] = ta[2];
    xx = r[0];
    yy = r[1];
    zz = r[2];
    switch ((int)i99632) {
	case 0: goto L14200;
    }
L15100:
    w = (rqsum - rho[1] * 2 * xx) / (sigma[1] * 2);
    denom = rr2 - xx * xx;
    sqr = rr2 - xx * xx - w * w;
    if (sqr < (float)0.) {
	if (dabs(sqr) > (float).001) {
	    s_wsfi(&io___309);
	    do_fio(&c__1, " Warning from CLSCHN -- Omega2 square is ", 41L);
	    do_fio(&c__1, (char *)&sqr, (ftnlen)sizeof(real));
	    e_wsfi();
	    cprint_(buffer, 200L);
	}
	sqr = (float)0.;
    }
    sqr = sqrt(sqr);
    comega[1] = (yy * w + sgns[0] * zz * sqr) / denom;
    somega[1] = (zz * w - sgns[0] * yy * sqr) / denom;
    q3[0] = xx - rho[1];
    q3[1] = w - sigma[1];
    q3[2] = sgns[0] * sqr;
    switch ((int)i99630) {
	case 0: goto L14300;
    }
L15200:
    comega[3] = (rho[2] * ctheta[3] - q3[0] * ctheta[2] - q3[1] * stheta[2]) /
	     (sigma[2] * stheta[3]);
/* Computing 2nd power */
    r__1 = comega[3];
    sqr = (float)1. - r__1 * r__1;
    if (sqr < (float)0.) {
	if (dabs(sqr) > (float).001) {
	    s_wsfi(&io___311);
	    do_fio(&c__1, " Warning from CLSCHN -- Omega4 square is ", 41L);
	    do_fio(&c__1, (char *)&sqr, (ftnlen)sizeof(real));
	    e_wsfi();
	    cprint_(buffer, 200L);
	}
	sqr = (float)0.;
	comega[3] = r_sign(&c_b366, &comega[3]);
    }
    sqr = sqrt(sqr);
    somega[3] = sgns[1] * sqr;
    switch ((int)i99628) {
	case 0: goto L14400;
    }
L15300:
    a = rho[2] * stheta[3] + sigma[2] * comega[3] * ctheta[3];
    b = q3[1] * ctheta[2] - q3[0] * stheta[2];
    c = sigma[2] * somega[3];
/* Computing 2nd power */
    r__1 = b;
/* Computing 2nd power */
    r__2 = q3[2];
    denom = r__1 * r__1 + r__2 * r__2;
    comega[2] = (a * b + q3[2] * c) / denom;
    somega[2] = (a * q3[2] - b * c) / denom;
    switch ((int)i99626) {
	case 0: goto L14500;
    }
L15400:
    ta[0] = ctheta[4];
    ta[1] = stheta[4];
    ta[2] = (float)0.;
    for (i = 4; i >= 1; --i) {
	tb[0] = ta[0];
	tb[1] = comega[i - 1] * ta[1] - somega[i - 1] * ta[2];
	tb[2] = somega[i - 1] * ta[1] + comega[i - 1] * ta[2];
	ta[0] = ctheta[i - 1] * tb[0] - stheta[i - 1] * tb[1];
	ta[1] = stheta[i - 1] * tb[0] + ctheta[i - 1] * tb[1];
	ta[2] = tb[2];
/* L15500: */
    }
    gphi = dot_(u, ta, &c__3) - ctheta[5];
    switch ((int)i99624) {
	case 0: goto L14600;
    }
L15600:
    if (omega[1] == lowphi || omega[1] == hiphi) {
	i99969 = 7;
	i99969_fmt = fmt_15700;
	goto L14100;
    }
L15700:
    omega[2] = atan2(somega[1], comega[1]);
    omega[3] = atan2(somega[2], comega[2]);
    omega[4] = atan2(somega[3], comega[3]);
    i99612 = 0;
    i99612_fmt = fmt_15800;
    goto L16300;
L15800:
    i99610 = 0;
    i99610_fmt = fmt_15900;
    goto L16500;
L15900:
    i99608 = 0;
    i99608_fmt = fmt_16000;
    goto L16800;
L16000:
    for (i = 1; i <= 5; i += 2) {
	if (pomega[(i + 1) / 2] != (float)0.) {
	    omega[i] += 3.141592653589794;
	}
/* L16100: */
    }
    for (i = 1; i <= 6; ++i) {
	d__1 = omega[i] / 6.283185307179588;
	omega[i] -= i_dnnt(&d__1) * 6.283185307179588;
/* L16200: */
    }
    switch ((int)i99856) {
	case 0: goto L4300;
    }
L16300:
    ta[0] = ctheta[0] * u[0] + stheta[0] * u[1];
    ta[1] = -(doublereal)stheta[0] * u[0] + ctheta[0] * u[1];
    ta[2] = u[2];
    for (i = 1; i <= 4; ++i) {
	i99602 = 0;
	i99602_fmt = fmt_16400;
	goto L16700;
L16400:
	;
    }
    if (debug > 0) {
	s_wsfi(&io___319);
	do_fio(&c__1, (char *)&ta[0], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ctheta[5], (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    comega[4] = ta[1] / stheta[5];
    somega[4] = ta[2] / stheta[5];
    omega[5] = atan2(somega[4], comega[4]);
    switch ((int)i99612) {
	case 0: goto L15800;
    }
L16500:
    ta[0] = ctheta[0] * v[0] + stheta[0] * v[1];
    ta[1] = -(doublereal)stheta[0] * v[0] + ctheta[0] * v[1];
    ta[2] = v[2];
    for (i = 1; i <= 5; ++i) {
	i99602 = 1;
	i99602_fmt = fmt_16600;
	goto L16700;
L16600:
	;
    }
    if (debug > 0) {
	s_wsfi(&io___320);
	do_fio(&c__1, (char *)&ta[0], (ftnlen)sizeof(real));
	e_wsfi();
	cprint_(buffer, 200L);
    }
    comega[5] = ta[1];
    somega[5] = ta[2];
    omega[6] = atan2(somega[5], comega[5]);
    switch ((int)i99610) {
	case 0: goto L15900;
    }
L16700:
    tb[0] = ta[0];
    tb[1] = comega[i - 1] * ta[1] + somega[i - 1] * ta[2];
    tb[2] = -(doublereal)somega[i - 1] * ta[1] + comega[i - 1] * ta[2];
    ta[0] = ctheta[i] * tb[0] + stheta[i] * tb[1];
    ta[1] = -(doublereal)stheta[i] * tb[0] + ctheta[i] * tb[1];
    ta[2] = tb[2];
    switch ((int)i99602) {
	case 0: goto L16400;
	case 1: goto L16600;
    }
L16800:
    for (i = 6; i >= 0; --i) {
	if (i != 0) {
	    newx[i] = (float)0.;
	    newy[i] = (float)0.;
	    newz[i] = (float)0.;
	}
	for (j = i + 1; j <= 6; ++j) {
	    ta[0] = newx[j];
	    ta[1] = newy[j] * comega[i] - newz[j] * somega[i];
	    ta[2] = newy[j] * somega[i] + newz[j] * comega[i];
	    newx[j] = ta[0] * ctheta[i] - ta[1] * stheta[i] + p[i * 3];
	    newy[j] = ta[0] * stheta[i] + ta[1] * ctheta[i] + p[i * 3 + 1];
	    newz[j] = ta[2] + p[i * 3 + 2];
/* L16850: */
	}
/* L16900: */
    }
    ind = atmind[4];
    for (i = 1; i <= 6; ++i) {
	ta[0] = x[ind];
	ta[1] = y[ind];
	ta[2] = z[ind];
	for (j = 1; j <= 3; ++j) {
	    ta[j - 1] = ta[j - 1] + newx[i] * u0[j - 1] + newy[i] * v0[j - 1] 
		    + newz[i] * w0[j - 1];
/* L16950: */
	}
	newx[i] = ta[0];
	newy[i] = ta[1];
	newz[i] = ta[2];
/* L17000: */
    }
    switch ((int)i99608) {
	case 0: goto L16000;
    }

/* OMUPD rkw 13/08/92 swapped inappropriate '/' with ',' */

/* 9001 FORMAT ('0Warning from CLSCHN -- Boundaries failed to match.'/ */
/*    +        ' ILIM = ',I1,' NP = ',I1/' G(PHI1) =',4(1X,1PG14.7)) */

} /* clschna_ */

