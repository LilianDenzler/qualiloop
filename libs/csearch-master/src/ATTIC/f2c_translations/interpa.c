/* interpa.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

doublereal interpa_(x1, y1, x2, y2, x)
real *x1, *y1, *x2, *y2, *x;
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static doublereal m;


/*     Performs linear interpolation between the points X1,Y1 and X2,Y2 */
/*     for X. The value of the INTERPOLATE is Y. The interpolation */
/*     formula attempts to postpone until the last moment the subtraction 
*/
/*     of roughly equal magnitude quantities. */

/*      IMPLICIT INTEGER(A-Z) */

    if (*x == *x1) {
	ret_val = *y1;
    } else if (*x == *y2) {
	ret_val = *y2;
    } else {
	m = ((doublereal) (*y2) - *y1) / ((doublereal) (*x2) - *x1);
	ret_val = m * ((doublereal) (*x) - *x1) + *y1;
    }
    return ret_val;
} /* interpa_ */

