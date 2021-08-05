/*****************************************************************************
 *      Name: adatmtgrid.c                                                   *
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
*   Adds atom ind to the space grid. If the coordinates of the atom
*   are out of bounds for the grid, then the grid is enlarged.
*/
void adatmtgrd(
int ind
)
{
   int ix,iy,iz,ispace,ispace1,p;
   int ind1;
   logical outofbound;
   int indgrid;
 
   if (grid.ingrid[ind1 = ind-1])
   {
      fprintf(out,"\nError in adatmtgrd -- \
Atom %d added twice to close contacts.\n",ind);
      die();
   }
   do
   {
      compspcgrd(&ind,&ix,&iy,&iz,&outofbound);
      if (outofbound) regrid(ind);
   } while (outofbound);
   indgrid = ((iz-1)*grid.ngridy+iy-1)*grid.ngridx+ix-1;
   ispace = grid.space_grid[indgrid];
   if (ispace == 0)
   {
      if (grid.freehd == 0)
      {
         fprintf(out,"\nError in adatmtgrd -- \
No more atom heads\n");
         die();
      }
      ispace = grid.freehd;
      ispace1 = ispace - 1;
      grid.freehd = grid.nexthd[grid.freehd-1];
      grid.space_grid[indgrid] = ispace;
      grid.clshd[ispace1] = 0;
      grid.clstl[ispace1] = 0;
   }
   else
      ispace1 = ispace - 1;
   if (grid.freecls == 0)
   {
      fprintf(out,"\nError in adatmtgrd -- No more close atoms available.\n");
      fprintf(out,"ind = %d\n",ind);
      dump_grid();
      die();
   }
   grid.clsatm[grid.freecls-1] = ind;
   p = grid.freecls;
   grid.freecls = grid.nextcls[grid.freecls-1];
   adafls(&grid.clshd[ispace1],&grid.clstl[ispace1],
          &grid.clstl[ispace1],&p,grid.nextcls);
   grid.ingrid[ind1] = true;
}
 
