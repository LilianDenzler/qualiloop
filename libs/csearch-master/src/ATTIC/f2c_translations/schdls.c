/* schdls.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"


/*     Copyright (c) 1987 Robert E. Bruccoleri */
/*     Copying of this software, in whole or in part, is permitted */
/*     provided that the copies are not made for commercial purposes, */
/*     appropriate credit for the use of the software is given, this */
/*     copyright notice appears, and notice is given that the copying */
/*     is by permission of Robert E. Bruccoleri. Any other copying */
/*     requires specific permission. */

/* Subroutine */ int schdls_(head, tail, free, next, a, item)
integer *head, *tail, *free, *next, *a, *item;
{
    extern /* Subroutine */ int trace_();
    static logical found;
    static integer itemp, prevp;
    extern /* Subroutine */ int cprint_();


/*     Searches for ITEM in the linked list specified by HEAD, TAIL, and 
*/
/*     NEXT whose elements are stored in A. The elements is deleted and */
/*     placed on the free list given by FREE. If the item is not found, */
/*     an error message is generated with a traceback, but execution will 
*/
/*     continue. */

/*      IMPLICIT INTEGER(A-Z) */
/* om omupd bnj 29/08/91 */
/* om */

/*     CRAP-OUT */
    /* Parameter adjustments */
    --a;
    --next;

    /* Function Body */
    if (*head != 0) {
	if (a[*head] == *item) {
	    itemp = *head;
	    if (*head == *tail) {
		*head = 0;
		*tail = 0;
	    } else {
		*head = next[*head];
	    }
	} else {
	    itemp = *head;
L20:
	    prevp = itemp;
	    itemp = next[itemp];
/*           CRAP-OUT */
	    if (itemp == 0) {
		goto L100;
	    }
	    found = a[itemp] == *item;
	    if (! found) {
		goto L20;
	    }
	    if (itemp == *tail) {
		*tail = prevp;
		next[*tail] = 0;
	    } else {
		next[prevp] = next[itemp];
	    }
	}
	next[itemp] = *free;
	*free = itemp;
	return 0;
    }
/* CC   TO CRAP-OUT */
L100:
    cprint_("***** Error in SCHDLS *****  Unable to find ITEM.", 49L);
    trace_();
    return 0;
} /* schdls_ */

