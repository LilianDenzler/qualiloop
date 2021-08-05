/* getparbond.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/*     Copyright (c) 1987 Robert E. Bruccoleri */
/*     Copying of this software, in whole or in part, is permitted */
/*     provided that the copies are not made for commercial purposes, */
/*     appropriate credit for the use of the software is given, this */
/*     copyright notice appears, and notice is given that the copying */
/*     is by permission of Robert E. Bruccoleri. Any other copying */
/*     requires specific permission. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal getparbond_(ib, jb, parm_no__, iac, kcb, ncb, type, cbb)
integer *ib, *jb;
shortint *parm_no__, *iac;
integer *kcb, *ncb, *type;
real *cbb;
{
    /* Format strings */
    static char fmt_9001[] = "(\002Error in GETPARBOND -- Unable to find par\
ameters for\002,\002 bond between \002,a4,\002 and \002,a4)";

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    static integer code, i, j;
    extern integer nindx_();
    static char buffer[100];
    extern /* Subroutine */ int cprint_(), die_();
    static integer icbb;

    /* Fortran I/O blocks */
    static icilist io___6 = { 0, buffer, 0, fmt_9001, 100, 1 };



/*     Returns the bond length for the bond between atoms IB and JB. */


    /* Parameter adjustments */
    --cbb;
    --type;
    --kcb;
    --iac;
    parm_no__ -= 101;

    /* Function Body */
    if (*ib == 0 || *jb == 0) {
	cprint_("Error in GETPARBOND: Atom index is zero", 39L);
	die_();
    }
    i = iac[*ib];
    j = iac[*jb];
    code = parm_no__[i + j * 100];
    icbb = nindx_(&code, &kcb[1], ncb);
    if (icbb == 0) {
	s_wsfi(&io___6);
	do_fio(&c__1, (char *)&type[*ib], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&type[*jb], (ftnlen)sizeof(integer));
	e_wsfi();
	cprint_(buffer, 100L);
	die_();
    }
    ret_val = cbb[icbb];
    return ret_val;
} /* getparbond_ */

