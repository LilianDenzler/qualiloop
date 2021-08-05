/* fortrace.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     FORTRACE */
/*     ACRM 21.06.91 */
/*     FORTRAN trace routine (VAX specific) */
/* Subroutine */ int fortrace_()
{
    extern /* Subroutine */ int cprint_();

    cprint_(" ", 1L);
    cprint_("A traceback will be generated via lib$signal()", 46L);
/*      CALL LIB$SIGNAL(%VAL(2)) */
    return 0;
} /* fortrace_ */

