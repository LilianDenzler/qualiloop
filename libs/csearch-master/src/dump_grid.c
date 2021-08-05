/*****************************************************************************
 *      Name: dump_grid.c                                                    *
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
 
/* Dump the contents of the space grid. */
 
void dump_grid(
void
)
{
   int ix,iy,iz,i,l,indgrid;
   char buf[500];
 
   for (iz=1; iz <= grid.ngridz; iz++)
   {
      fprintf(out,"\nSpace grid z = %d\n",iz);
      for (iy=1; iy <= grid.ngridy; iy++)
      {
         buf[0] = '\0';
         l = 0;
         for (ix=1; ix <= grid.ngridx; ix++)
         {
            indgrid = ((iz-1)*grid.ngridy+iy-1)*grid.ngridx+ix-1;
            l += sprintf(&buf[l],"%5d",grid.space_grid[indgrid]);
         }
         fprintf(out,"%4d: %s\n",iy,buf);
      }
   }
   fprintf(out,"\nfreecls = %d  freehd = %d\n",grid.freecls,grid.freehd);
   fprintf(out,"\n    i    ingrid   nextcls    clsatm     clshd     clstl    \
nexthd\n");
   for (i=1; i<= values.natoms; i++)
      fprintf(out,"%5d%10d%10d%10d%10d%10d%10d\n",
             i,grid.ingrid[i-1],grid.nextcls[i-1],grid.clsatm[i-1],
             grid.clshd[i-1],grid.clstl[i-1],grid.nexthd[i-1]);
}
 
