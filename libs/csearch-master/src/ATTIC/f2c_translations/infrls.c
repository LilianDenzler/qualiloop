/* infrls.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int infrls_(head, tail, next, n)
integer *head, *tail, *next, *n;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;


/*     Initializes an array into a free list, i.e. everything in the */
/*     array is linked sequentially into one long list. The length of the 
*/
/*     array is given by N. */

/*      IMPLICIT INTEGER(A-Z) */
/* OM OMUPD BNJ 29/08/91 */
/* OM */

    /* Parameter adjustments */
    --next;

    /* Function Body */
    if (*n == 0) {
	*head = 0;
	*tail = 0;
    } else {
	*head = 1;
	*tail = *n;
	i__1 = *n;
	for (i = 2; i <= i__1; ++i) {
	    next[i - 1] = i;
/* L50: */
	}
	next[*n] = 0;
    }
    return 0;
} /* infrls_ */

