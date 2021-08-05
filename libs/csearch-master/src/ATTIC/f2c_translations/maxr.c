/* maxr.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     Recoded ACRM 12.06.91 */
doublereal maxr_(array, num)
real *array;
integer *num;
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i;

    /* Parameter adjustments */
    --array;

    /* Function Body */
    ret_val = array[1];
    i__1 = *num;
    for (i = 2; i <= i__1; ++i) {
	if (array[i] > ret_val) {
	    ret_val = array[i];
	}
/* L100: */
    }
    return ret_val;
} /* maxr_ */

