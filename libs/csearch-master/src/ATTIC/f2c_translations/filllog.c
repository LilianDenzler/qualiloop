/* filllog.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     FILLLOG */
/*     ACRM 24.06.91 */
/*     Fills an array with logicals */
/* Subroutine */ int filllog_(arr, num, truth)
logical *arr;
integer *num;
logical *truth;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;

    /* Parameter adjustments */
    --arr;

    /* Function Body */
    i__1 = *num;
    for (i = 1; i <= i__1; ++i) {
	arr[i] = *truth;
/* L100: */
    }
    return 0;
} /* filllog_ */

