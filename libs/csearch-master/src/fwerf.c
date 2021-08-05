/*****************************************************************************
 *      Name: fwerf                                                          *
 *  Function: to call system erf routine and typecast output                 *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1992                                 *
 *---------------------------------------------------------------------------*
 *    Author: Brian Jamieson                                                 *
 *      Date: 16/03/92                                                       *
 *---------------------------------------------------------------------------*
 *    Inputs: *xx, *xtemp                                                    *
 *   Outputs:                                                                *
 *   Returns:                                                                *
 * Externals:                                                                *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/

#include "ProtoTypes.h"
#include "CongenProto.h"
extern double erf(double);

void fwerf(float *xx, float *xtemp)
{
double dxx,dxtemp;

   dxx = (double) (*xx);

   dxtemp = erf(dxx);
   /* printf("erf of %f is %f\n",dxx,dxtemp); */

   *xtemp = (float) dxtemp;
}
