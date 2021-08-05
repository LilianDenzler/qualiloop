/* cprint.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     CPRINT */
/*     ACRM 13.06.91 */
/*     Replaces writes to channel 6 with C i/o */
/* Subroutine */ int cprint_(string, string_len)
char *string;
ftnlen string_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len();

    /* Local variables */
    extern /* Subroutine */ int cprint1_();

/*     Modified for f2c    ACRM 07.01.94 */
/*      CALL CPRINT1(%REF(STRING),LEN(STRING)) */
    i__1 = i_len(string, string_len);
    cprint1_(string, &i__1);
    return 0;
} /* cprint_ */

