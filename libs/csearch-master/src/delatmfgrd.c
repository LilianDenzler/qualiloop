/*****************************************************************************
 *      Name: delatmfgrd.c                                                   *
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
 
void delatmfgrd(
int ind
)
{
   delatmfgrd1(&ind,grid.space_grid,grid.ingrid,
                          grid.nexthd,grid.clshd,grid.clstl,
                          grid.nextcls,grid.clsatm);
}
 
