/*****************************************************************************
 *      Name: extnlims.c                                                     *
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
#include "values.h"
 
void extnlims(
void
)
{
   grid.xmx = maxcoor(coords.xcart,values.natoms);
   grid.ymx = maxcoor(coords.ycart,values.natoms);
   grid.zmx = maxcoor(coords.zcart,values.natoms);
   grid.xmn = mincoor(coords.xcart,values.natoms);
   grid.ymn = mincoor(coords.ycart,values.natoms);
   grid.zmn = mincoor(coords.zcart,values.natoms);
   prtgrdsz();
}
 
