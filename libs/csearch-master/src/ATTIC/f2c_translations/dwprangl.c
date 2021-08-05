/* dwprangl.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int dwprangl_(it, jt, kt, parm_no__, iac, kct, nct, natc, 
	type, ctb, retval)
integer *it, *jt, *kt;
shortint *parm_no__, *iac;
integer *kct, *nct, *natc, *type;
real *ctb, *retval;
{
    /* Format strings */
    static char fmt_9001[] = "(\002Error in GTPRANGL: Atom index is zero\002)"
	    ;
    static char fmt_9002[] = "(\0020Error in GTPRANGL -- Unable to find para\
meters for\002,\002 angle between \002,a4,\002 and \002,a4,\002 and \002,a4)";

    /* Builtin functions */
    integer s_wsfi(), e_wsfi(), do_fio();

    /* Local variables */
    static integer code, ictb, i, j, k;
    extern integer nindx_();
    static char buffer[100];
    extern /* Subroutine */ int cprint_(), die_();

    /* Fortran I/O blocks */
    static icilist io___2 = { 0, buffer, 0, fmt_9001, 100, 1 };
    static icilist io___8 = { 0, buffer, 0, fmt_9002, 100, 1 };



/*     Returns the bond angle for the angle between atoms IT, JT, & KT. */


    /* Parameter adjustments */
    --ctb;
    --type;
    --kct;
    --iac;
    parm_no__ -= 101;

    /* Function Body */
    if (*it == 0 || *jt == 0 || *kt == 0) {
	s_wsfi(&io___2);
	e_wsfi();
	cprint_(buffer, 100L);
	die_();
    }
    i = iac[*it];
    j = iac[*jt];
    k = iac[*kt];
    code = parm_no__[i + k * 100];
    code = *natc * code + j;
    ictb = nindx_(&code, &kct[1], nct);
    if (ictb == 0) {
	s_wsfi(&io___8);
	do_fio(&c__1, (char *)&type[*it], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&type[*jt], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&type[*kt], (ftnlen)sizeof(integer));
	e_wsfi();
	cprint_(buffer, 100L);
	die_();
    }
    *retval = ctb[ictb];
    return 0;
} /* dwprangl_ */

