/*****************************************************************************
 *      Name: fill_grid.c                                                    *
 *  Function:                                                                *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1992                                 *
 *---------------------------------------------------------------------------*
 *    Author: B. N. Jamieson, Oxford Molecular  Ltd                          *
 *      Date: 26/02/92                                                       *
 *---------------------------------------------------------------------------*
 *    Inputs:                                                                *
 *   Outputs:                                                                *
 *   Returns:                                                                *
 * Externals:                                                                *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 * 26/02/92   DW         Added header; removed invalid "grid.H"              *
 *****************************************************************************/
#include "ProtoTypes.h"
#include "CongenProto.h"
 
void fill_grid(
void
)
{
   fill_grid1(grid.space_grid,grid.excluded,grid.cntnbx,grid.ingrid,
              grid.nbxa,grid.nexthd,grid.clshd,grid.clstl,grid.nextcls,
              grid.clsatm,grid.resbya,grid.radius);
}
 
