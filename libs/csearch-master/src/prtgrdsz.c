/*****************************************************************************
 *      Name: prtgrdsz.c                                                     *
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
 
/* Prints the current size of the space grid. */
 
void prtgrdsz(
void
)
{
   fprintf(out,"\nLimits of box around molecule:\n");
   fprintf(out,"\nX: %10.5f  %10.5f\n",grid.xmn,grid.xmx);
   fprintf(out,"Y: %10.5f  %10.5f\n",grid.ymn,grid.ymx);
   fprintf(out,"Z: %10.5f  %10.5f\n",grid.zmn,grid.zmx);
}
 
