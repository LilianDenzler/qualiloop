/*****************************************************************************
 *      Name: getgrdspc.c                                                    *
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
 
void getgrdspc (
void
)
{
   float f;
 
   grid.spgridsz = 2.75;
   grid.recipgrid = 1.0 / grid.spgridsz;
   /* ACRM added casts */
   grid.ngridx = chmceil((f=(float)((grid.xmx-grid.xmn)*grid.recipgrid),&f));
   grid.ngridy = chmceil((f=(float)((grid.ymx-grid.ymn)*grid.recipgrid),&f));
   grid.ngridz = chmceil((f=(float)((grid.zmx-grid.zmn)*grid.recipgrid),&f));
 
   grid.space_grid = (short int *)
             alloc(grid.ngridx*grid.ngridy*grid.ngridz*sizeof(short int));
   grid.maxnbx = 15;
   grid.clshd = (int *) alloc(values.natoms*sizeof(int));
   grid.clstl = (int *) alloc(values.natoms*sizeof(int));
   grid.clsatm = (int *) alloc(values.natoms*sizeof(int));
   grid.nextcls = (int *) alloc(values.natoms*sizeof(int));
   grid.nexthd = (int *) alloc(values.natoms*sizeof(int));
   grid.nbxa = (short int *) alloc(grid.maxnbx*values.natoms*sizeof(short int));
   grid.cntnbx = (int *) alloc(values.natoms*sizeof(int));
   grid.excluded = (logical *) alloc(values.natoms*sizeof(logical));
   grid.donp = (int *) alloc(values.natoms*sizeof(int));
   grid.accp = (int *) alloc(values.natoms*sizeof(int));
   grid.ingrid = (logical *) alloc(values.natoms*sizeof(logical));
   grid.resbya = (int *) alloc(values.natoms * sizeof(int));
   /* ACRM Corrected sizeof(int) to sizeof(logical) */
   grid.qside = (logical *) alloc(values.natoms * sizeof(logical));
   grid.radius = (float *) alloc(values.natoms * sizeof(float));
}
 
