/* abmerf.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/****************************************************************************
***/
/**      Name: abmerf                                                         
**/
/**  Function: to wrap and use the system erf function                        
**/
/** Copyright: (C) OXFORD MOLECULAR LTD, 1992                                 
**/
/**--------------------------------------------------------------------------
-**/
/**    Author: Brian Jamieson                                                 
**/
/**      Date: 16/03/92                                                       
**/
/**--------------------------------------------------------------------------
-**/
/**    Inputs: x                                                             
**/
/**   Outputs:                                                                
**/
/**   Returns: abmerf                                                         
**/
/** Externals:                                                                
**/
/**--------------------------------------------------------------------------
-**/
/** MODIFICATION RECORD:                                                      
**/
/** DD/MM/YY   Initials   Comments                                            
**/
/****************************************************************************
*+/*/

doublereal abmerf_(x)
real *x;
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    extern /* Subroutine */ int fwerf_();
    static real xtemp, xx;




    xx = *x;

    fwerf_(&xx, &xtemp);

    ret_val = xtemp;

    return ret_val;
} /* abmerf_ */

