/* nindx.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     NINDX */
/*     Recoded ACRM 12.06.91 */
/*     Finds NUMBER in sorted NARRAY (length NLEN) by binary search */
/*     Returns its index in the array, 0 if not found. */
integer nindx_(number, narray, nlen)
integer *number, *narray, *nlen;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer istop, istart;

    /* Parameter adjustments */
    --narray;

    /* Function Body */
    istart = 1;
    istop = *nlen;
L100:
    ret_val = istart + (istop - istart) / 2;
    if (narray[ret_val] == *number) {
	return ret_val;
    } else if (narray[ret_val] > *number) {
	istop = ret_val - 1;
    } else if (narray[ret_val] < *number) {
	istart = ret_val + 1;
    }
    if (istop >= istart) {
	goto L100;
    }
    ret_val = 0;
    return ret_val;
} /* nindx_ */

