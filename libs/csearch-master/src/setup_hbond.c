/*****************************************************************************
 *      Name: setup_hbond.c                                                  *
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
*   Sets up various variables and arrays for computing the hydrogen
*   bond energy in sidechain atoms, and also takes care of HB variables
*   stored with each sidechain degree of freedom.
*/
 
void setup_hbond(
void
)
{
   int i,*ip,*jp,n;
   struct dof *p;
   struct sidechain_d *desc;
   struct sideres *srp,**srpp;
 
   stuphb1(grid.donp,grid.accp);
   for (p=dof_head; p; p = p->next)
   {
      if (p->dof_type == sidechain_t)
      {
         desc = (struct sidechain_d *)p->desc; /* ACRM added cast */
         for (srpp=desc->residues; (srp = *srpp); srpp++)
         {
            n = 0;
            for (ip=srp->atomp; (i = *ip); ip++)
               if (grid.donp[i-1]) n++;
            jp = srp->donors = (int *) alloc((n+1)*sizeof(int));
            for (ip=srp->atomp; (i = *ip); ip++)
               if (grid.donp[i-1]) *jp++ = i;
            *jp = 0;
            n = 0;
            for (ip=srp->atomp; (i = *ip); ip++)
               if (grid.accp[i-1]) n++;
            jp = srp->acceptors = (int *) alloc((n+1)*sizeof(int));
            for (ip=srp->atomp; (i = *ip); ip++)
               if (grid.accp[i-1]) *jp++ = i;
            *jp = 0;
         }
      }
   }
}
 
