/* forprtatm.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     FORPRTATM */
/*     ACRM 18.08.91 */
/*     Fortran call to print_atom1() */
/* Subroutine */ int forprtatm_(string, atnum, string_len)
char *string;
integer *atnum;
ftnlen string_len;
{
    /* Builtin functions */
    /* Subroutine */ int s_copy();

    /* Local variables */
    extern /* Subroutine */ int c_print_atom__();

    s_copy(string, " ", string_len, 1L);
/* Modified for f2c */
/*      CALL C_PRINT_ATOM(%REF(STRING),%VAL(ATNUM)) */
    c_print_atom__(string, atnum);
    return 0;
} /* forprtatm_ */

