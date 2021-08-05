/*****************************************************************************
 *      Name: schnatm.c                                                      *
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
*   In the event atom ind doesn't fit in the grid, the grid is enlarged.
*/
 
void schnatm(
int ind,
int mode,
float *maxevdw,
logical ignore_evdw,
int *sidehits,
logical *impact,
float *eel,
float *evdw,
float *ehb
)
{
   logical outofbound;
/* ACRM added these for call to SRCHNAT1 */
   double eel_d,
          evdw_d,
          ehb_d;
 
   if (mode == contact && *maxevdw >= 1.0e20)
      *impact = false;
   else
   {
      do
      {
/* ACRM Corrected to doubles
//         schnatm1(&ind,grid.space_grid,&mode,maxevdw,&ignore_evdw,
//                           grid.cntnbx,grid.nbxa,
//                           grid.excluded,cg.parm_no,grid.clshd,grid.nextcls,
//                           grid.clsatm,grid.donp,grid.accp,grid.qside,sidehits,
//                           grid.resbya,impact,eel,evdw,ehb,&outofbound);
*/
         eel_d  = *eel;
         evdw_d = *evdw;
         ehb_d  = *ehb;
         srchnatm1(&ind,grid.space_grid,&mode,maxevdw,&ignore_evdw,
                           grid.cntnbx,grid.nbxa,
                           grid.excluded,cg.parm_no,grid.clshd,grid.nextcls,
                           grid.clsatm,grid.donp,grid.accp,grid.qside,sidehits,
                           grid.resbya,impact,&eel_d,&evdw_d,&ehb_d,
                           &outofbound);
         *eel  = eel_d;
         *evdw = evdw_d;
         *ehb  = ehb_d;
 
         if (outofbound) regrid(ind);
      }  while (outofbound);
   }
}
 
