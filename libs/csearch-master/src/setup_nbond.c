/*****************************************************************************
 *      Name: setup_nbond.c                                                  *
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
 
/* Sets up data needed for computing the non-bonded interaction. */
 
void setup_nbond(
void
)
{
   int itri,i;
   struct dof *dofp;
   struct atom **atpp;
 
/* ACRM
 
   CONGEN stores all nbond info on the heap. bases.bnbnd is an
   array of heap addresses used to store the information while
   engpar.* is a set of constant parameter indices into bases.bnbnd
   which tell you where to look for a value.
   We will change this so that engpar.* is where stuff is actually
   stored.
*/
 
/*
// cg.cutnb2 = rheap[bases.bnbnd[engpar.nbcut-1]-1];
// cg.cutnb2 *= cg.cutnb2;
*/
   cg.cutnb2 = engpar.nbcut * engpar.nbcut;
 
 
/*
// cg.epsilon = rheap[bases.bnbnd[engpar.dielec-1]-1];
*/
   cg.epsilon = engpar.dielec;
 
/* Constant dielectric always used
// cg.cons_die = logical_f(heap[bases.bnbnd[engpar.nbflag-1]-1] == 7);
*/
   cg.cons_die = f77_true;
 
   for (itri=0,i=1; i<=100; itri += i,i++) cg.ioff[i-1] = itri;
   cg.grid2 = grid.spgridsz * grid.spgridsz;
   /* ACRM This use of fill4 is with 'logicals' (i.e. int's)
      so should be OK providing int is same as long */
   fill4((long *)grid.qside,&values.natoms,(long *)(i=f77_false,&i));
   for (dofp = dof_head; dofp; dofp = dofp->next)
      if (dofp->dof_type == sidechain_t)
         for (atpp = dofp->atoms; *atpp; atpp++)
            grid.qside[(*atpp)->atomno-1] = f77_true;
}
 
