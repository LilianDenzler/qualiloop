/*****************************************************************************
 *      Name: getclsangle.c                                                  *
 *  Function:                                                                *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1992                                 *
 *---------------------------------------------------------------------------*
 *    Author: B.N. Jamieson, Oxford Molecular Ltd                            *
 *      Date: 26/02/92                                                       *
 *---------------------------------------------------------------------------*
 *    Inputs:                                                                *
 *   Outputs:                                                                *
 *   Returns:                                                                *
 * Externals:                                                                *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 * 26/02/92   DW         Added header; removed grid.H                        *
 *****************************************************************************/

#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Gets the bond angle between atoms i, j, and k as required by chain
*   closure. If any of the atoms in the angle are constructed in this
*   run of CGEN, then the parameters are used for the angle.
*   Otherwise, the actual bond angle as given by the coordinates is
*   used.
*/
 
float getclsangle(
int i,
int j,
int k
)
{
   logical icons,jcons,kcons;
   short int si,sj,sk;
   float a;
   logical _false_ = f77_false;
 
   icons = !ref(grid.ingrid,i-1);
   jcons = !ref(grid.ingrid,j-1);
   kcons = !ref(grid.ingrid,k-1);
   if (icons || jcons || kcons)
      a = gtprangl2(i,j,k);
   else
   {
      /* ACRM added casts */
      angle((si=(short)i,&si),(sj=(short)j,&sj),(sk=(short)k,&sk),
            &_false_,&_false_,&one,&a,coords.xcart,coords.ycart,coords.zcart);
      a *= dtorad;
   }
   return a;
}
 
