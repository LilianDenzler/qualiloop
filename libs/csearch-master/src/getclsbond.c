/*****************************************************************************
 *      Name: getclsbond.c                                                   *
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
*   Gets the bond length between atoms i and j as required by chain closure.
*   If any of the atoms in the bond are constructed in this run of CGEN,
*   then the parameters are used for the bond length. Otherwise, the
*   actual bond length as given by the coordinates are used.
*/
 
float getclsbond(
int i,
int j
)
{
   logical icons,jcons;
   short int si,sj,sk,sone;
   float bond;
   /* ACRM corrected int to logical */
   logical _false_ = f77_false;
 
   /* ACRM ref() is just a macro which gives the same
      effect as grid.ingrid[i-1] or grid.ingrid[j-1]
   */
   icons = !ref(grid.ingrid,i-1);
   jcons = !ref(grid.ingrid,j-1);
   if (icons || jcons)
      bond = getparbond2(i,j);
   else
/*    ACRM Corrected to short for sk and sone; also added casts
//      bondl((si=i,&si),(sj=j,&sj),&zero,&_false_,&one,
//            coords.xcart,coords.ycart,coords.zcart,&bond);
*/
      bondl((si=(short)i,&si),(sj=(short)j,&sj),(sk=(short)zero,&sk),
            &_false_,(sone=(short)one,&sone),coords.xcart,coords.ycart,coords.zcart,&bond);
 
   return bond;
}
 
