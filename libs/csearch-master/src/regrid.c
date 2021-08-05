/*****************************************************************************
 *      Name: regrid.c                                                       *
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
 
/*
*   Rebuilds the space grid. In order to avoid any problems with
*   roundoff error, we rebuild the entire array rather than attempting
*   the moving of elements between the old and new arrays. Since
*   this rebuilding operation should not occur very often, this
*   strategy will not be very expensive.
*/
 
void regrid(
int ind
)
{
   int ind1;
   int i,dummy_tail,count;
   float f;
   short int s1;
 
   if (dbg.cgen)
   {
      fprintf(out,"Regridding: Old bounds\n");
      prtgrdsz();
   }
   ind1 = ind - 1;
   for (i=0; i<values.natoms; i++)
   {
      if (coords.xcart[i] != anum && (i == ind1 || grid.ingrid[i]) )
      {
         grid.xmn = minf2(grid.xmn,coords.xcart[i]);
         grid.ymn = minf2(grid.ymn,coords.ycart[i]);
         grid.zmn = minf2(grid.zmn,coords.zcart[i]);
         grid.xmx = maxf2(grid.xmx,coords.xcart[i]);
         grid.ymx = maxf2(grid.ymx,coords.ycart[i]);
         grid.zmx = maxf2(grid.zmx,coords.zcart[i]);
      }
   }
   /* Added casts */
   grid.ngridx = chmceil((f=(float)((grid.xmx-grid.xmn)*grid.recipgrid),&f));
   grid.ngridy = chmceil((f=(float)((grid.ymx-grid.ymn)*grid.recipgrid),&f));
   grid.ngridz = chmceil((f=(float)((grid.zmx-grid.zmn)*grid.recipgrid),&f));
   free(grid.space_grid);
   grid.space_grid = (short int *)
      alloc(grid.ngridx*grid.ngridy*grid.ngridz*sizeof(short int));
   /* Should be no problems here; added casts */
   fill2(grid.space_grid,(i=(int)(grid.ngridx*grid.ngridy*grid.ngridz),&i),
         (s1=(short)0,&s1));
   infrls(&grid.freehd,&dummy_tail,grid.nexthd,&values.natoms);
   infrls(&grid.freecls,&dummy_tail,grid.nextcls,&values.natoms);
   count = 0;
   for (i=1; i <= values.natoms; i++)
   {
      if (grid.ingrid[i-1])
      {
         grid.ingrid[i-1] = false;
         adatmtgrd(i);
         if (++count != grid.freecls-1)
         {
            fprintf(out,"Discrepancy noted in regrid.\n");
            fprintf(out,"Count = %d  i = %d  freecls = %d.\n",
                   count,i,grid.freecls);
         }
      }
   }
   if (dbg.cgen)
   {
      fprintf(out,"Regridding: New bounds\n");
      prtgrdsz();
   }
}
 
