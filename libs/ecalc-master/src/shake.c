/*************************************************************************

   Program:    ECalc
   File:       shake.c
   
   Version:    V1.5.1
   Date:       07.01.21
   Function:   Perform the SHAKE algorithm to regularize bond lengths
   
   Copyright:  (c) UCL, Prof. Andrew C. R. Martin 1995-2021
   Author:     Prof. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.4   18.05.95 Original
   V1.5   06.02.03 Skipped
   V1.5.1 07.01.21 Removed unused variables

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/macros.h"

#include "topology.h"
#include "ecalc.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/


/************************************************************************/
/*>int ShakeChain(TOPOLOGY *topol, REAL Tolerence)
   -----------------------------------------------
   Input:   TOPOLOGY  *topol        Topology pointer
            REAL      Tolerence     Required geometric tolerence
   Returns: int                     Number of iterations performed

   Performs the SHAKE algorithm to regularise bond lengths.

   18.05.95 Original    By: ACRM
   19.05.95 Errors reported through StoreError()
*/
int ShakeChain(TOPOLOGY *topol, REAL Tolerence)
{
   BOOL  Ready;
   BOND  *p;
   ATOM  *atom1, 
         *atom2;
   int   NIter,
         i;
   REAL  Tol2,
         OptLenSq,
         Diff,
         DistIJSq,
         DotProduct,
         InvMassI,
         InvMassJ,
         ResetDist,
         xh;
   VEC3F VectorIJ,
         VectorRefIJ;


   if(topol->NBonds == 0 || gNumBondParams == 0)
      return(1);
   
   /* Initialise variables                                              */
   NIter = 0;
   Tol2  = 2.E0*Tolerence;
   
   /* For each atom in the molecule set two flags                       */
   for(i=0; i<topol->NAtoms; i++)
   {
      (topol->atoms[i])->Skip     = TRUE;
      (topol->atoms[i])->SkipPrev = FALSE;
   }
   
   /* Clear flag which indicates completion of the loop                 */
   Ready=FALSE;
   
   /* Begin the loop                                                    */
   for(;;)
   {
      /* Break out of loop if completed                                 */
      if(Ready)
         return(NIter);
      
      /* Check iteration count                                          */
      if(NIter > MAXSHAKEITER) 
      {
         sprintf(gError,"SHAKE exceeded %d iterations\n",MAXSHAKEITER);
         StoreError("ShakeChain()",gError);
         
         return(0);
      }
      
      /* Assume that this will be our last iteration                    */
      Ready=TRUE;
      
      /* For each bond being shaken                                     */
      for(p=topol->bonds; p!=NULL; NEXT(p))
      {
         atom1 = p->atom1;
         atom2 = p->atom2;
         
         if(VALID(atom1) && VALID(atom2))
         {
            if(atom1->SkipPrev && atom2->SkipPrev)
               break;
            
            OptLenSq = p->OptLen * p->OptLen;
         }
         
         Diff     = OptLenSq;            
         
         /* Record the vector between the 2 atoms from the coord set
            to be shaken in VectorIJ and calculate the squared distance
            between then in DistIJSq 
         */
         VectorIJ.x = atom1->x - atom2->x;
         VectorIJ.y = atom1->y - atom2->y;
         VectorIJ.z = atom1->z - atom2->z;
         DistIJSq   = DISTSQ(atom1, atom2);
         
         /* See if we're within tolerence                               */
         Diff -= DistIJSq;
         
         if(ABS(Diff) < OptLenSq*Tol2)
            break;
         
         /* Record the vector between the 2 atoms from the reference
            set in VectorRefIJ and calculate the squared distance 
            between then in RIJ2 
         */
         VectorRefIJ.x = atom1->Refx - atom2->Refx;
         VectorRefIJ.y = atom1->Refy - atom2->Refy;
         VectorRefIJ.z = atom1->Refz - atom2->Refz;
         
         /* Calculate the dot product between the vectors of previous
            and current steps 
         */
         DotProduct = (VectorIJ.x * VectorRefIJ.x) +
                      (VectorIJ.y * VectorRefIJ.y) +
                      (VectorIJ.z * VectorRefIJ.z);

         /* If the dot product is too small print error message and 
            die 
         */
         if(DotProduct < OptLenSq*1.E-6) 
         {
            sprintf(gError, "SHAKE failed, deviation too large\n\
NIter: %d, Atom1: %s, Atom2: %s\n",
                    NIter,atom1->atom,atom2->atom);
            StoreError("ShakeChain()",gError);
            
            return(0);
         }
         
         /* Modify the coordinates based on the mass of the atoms       */
         InvMassI = (REAL)1.0/atom1->mass;
         InvMassJ = (REAL)1.0/atom2->mass;
         
         ResetDist=Diff/(DotProduct*(InvMassI+InvMassJ)*(REAL)2.0);
         
         xh = VectorRefIJ.x * ResetDist;
         atom1->x += xh * InvMassI;
         atom2->x -= xh * InvMassJ;
         
         xh = VectorRefIJ.y * ResetDist;
         atom1->y += xh * InvMassI;
         atom2->y -= xh * InvMassJ;
         
         xh = VectorRefIJ.z * ResetDist;
         atom1->z += xh * InvMassI;
         atom2->z -= xh * InvMassJ;
         
         /* Set flags to say that these atoms have been modified        */
         atom1->Skip = FALSE;            
         atom2->Skip = FALSE;            
         Ready       = FALSE;            
         
      }  /* Continue with next bond                                     */
      
      /* Increment iteration count                                      */
      NIter++;
      
      /* Copy current Skip list and set current flags to TRUE           */
      for(i=0; i<topol->NAtoms; i++)
      {
         (topol->atoms[i])->SkipPrev = (topol->atoms[i])->Skip;
         (topol->atoms[i])->Skip     = TRUE;
      }
   }
}


/************************************************************************/
/*>void InitRelax(MOLECULE *mol)
   -----------------------------
   Initialises for performing the relaxation by copying the coordinates
   into the reference set.

   19.05.95 Original    By: ACRM
*/
void InitRelax(MOLECULE *mol)
{
   int i,
       j;
   TOPOLOGY *topol;
   

   for(j=0; j<mol->NChains; j++)
   {
      topol = mol->topol[j];

      /* For each atom in the molecule set two flags                    */
      for(i=0; i<topol->NAtoms; i++)
      {
         (topol->atoms[i])->Refx = (topol->atoms[i])->x;
         (topol->atoms[i])->Refy = (topol->atoms[i])->y;
         (topol->atoms[i])->Refz = (topol->atoms[i])->z;
      }
   }
}

   
/************************************************************************/
/*>int ShakeAll(MOLECULE *mol, REAL ShakeTol)
   ------------------------------------------
   Input:   MOLECULE *mol        Molecule pointer
            REAL     ShakeTol    Shake tolerence
   Returns: int                  Mean number of iterations

   Applies the SHAKE algorithm to all chains in the topology.

   N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.

   N.B. Disulphides are not correctly treated at present as they exist
   as separate topologies. These need to be patched into the main 
   topology of each chain to work correctly.

   N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.N.B.

   19.05.95 Original    By: ACRM
*/
int ShakeAll(MOLECULE *mol, REAL ShakeTol)
{
   int chain,
       NChains = 0,
       NIter   = 0,
       NI;
   TOPOLOGY *topol;
   
   for(chain=0; chain<mol->NChains; chain++)
   {
      topol = mol->topol[chain];
      
      if(!topol->Disulphide)
      {
         if((NI = ShakeChain(topol, ShakeTol)) != 0)
         {
            NIter += NI;
            NChains++;
         }
      }
   }
   if(NChains)
      NIter = (int)(NIter / NChains);

   return(NIter);
}


/************************************************************************/
/*>void RelaxStructure(MOLECULE *mol, REAL VDWTol, REAL ShakeTol)
   --------------------------------------------------------------
   19.05.95 Original    By: ACRM
*/
void RelaxStructure(MOLECULE *mol, REAL VDWTol, REAL ShakeTol)
{
   int  NIter,
        NVDW;
        /* NShake;                                                      */
   
   
   /* Initialise by copying the coordinates into the reference arrays   */
   InitRelax(mol);

   for(NIter=0; NIter<MAXRELAXITER; NIter++)
   {
      NVDW   = VDWRelax(mol, VDWTol);
      /* NShake = */
      ShakeAll(mol, ShakeTol);
      
      if(!NVDW)
         break;

      InitRelax(mol);
   }
}


/************************************************************************/
/*>int VDWRelax(MOLECULE *mol, REAL VDWTol)
   ----------------------------------------
   Performs a VDW relaxation by moving atoms apart along the vector
   between their centres. This is done in an iterative fashion weighting
   the atoms' movements by their respective masses in an analagous manner
   to SHAKE.

   19.05.95 Original    By: ACRM
*/
int VDWRelax(MOLECULE *mol, REAL VDWTol)
{
   int   i,
         j,
         NIter;
   ATOM  *atom1,
         *atom2;
   REAL  DistSq,
         OptDist,
         OptDistSq,
         TotMass,
         Factor,
         DotProduct,
         sep,
         SepSq;
   BOOL  Finished;
   VEC3F Vector,
         MoveVector;


   if(gNumNonBondParams == 0) 
      return(0);

   if(gGrid == NULL)
   {
      StoreError("VDWRelax()","No non-bond grid has been defined");
      return(0);
   }

   for(NIter=0; NIter<MAXVDWITER; NIter++)
   {
      Finished = TRUE;
      
      for(i=0; i<mol->NAtoms; i++)
      {
         atom1 = gGrid[i].atom;
         if(VALID(atom1))
         {
            for(j=0; j<gGrid[i].NNeighb; j++)
            {
               atom2 = gGrid[i].neighbours[j];
               if(VALID(atom2))
               {
                  /* Find the optimum distance between the 2 atoms      */
                  OptDistSq = OptDist = atom1->vdwr + atom2->vdwr;
                  OptDistSq *= OptDistSq;
                  
                  /* If distance is less than the optimum distance, but is
                     less than the optimum distance - tolerence
                  */
                  if(((DistSq=DISTSQ(atom1, atom2)) != (REAL)0.0) &&
                     ((OptDistSq - DistSq) > (OptDistSq * VDWTol)))
                  {
                     if(atom1->Mobile || atom2->Mobile)
                     {
                        Finished = FALSE;
                     
                        /* Calculate the amount by which the atoms need 
                           to be separated
                        */
                        sep   = (REAL)0.1 *  
                                (OptDist - (REAL)sqrt((double)DistSq));
                        SepSq = sep * sep;

                        /* Calculate vector from atom1 to atom2         */
                        Vector.x = atom2->x - atom1->x;
                        Vector.y = atom2->y - atom1->y;
                        Vector.z = atom2->z - atom1->z;

                        DotProduct = (Vector.x * Vector.x) +
                                     (Vector.y * Vector.y) +
                                     (Vector.z * Vector.z);

                        Factor = SepSq / DotProduct;

                        MoveVector.x = sqrt(Factor * Vector.x * Vector.x);
                        MoveVector.y = sqrt(Factor * Vector.y * Vector.y);
                        MoveVector.z = sqrt(Factor * Vector.z * Vector.z);

                        if(atom1->Mobile && atom2->Mobile)
                        {
                           TotMass = atom1->mass + atom2->mass;
                        
                           atom1->x -= atom2->mass * MoveVector.x /
                                       TotMass;
                           atom1->y -= atom2->mass * MoveVector.y /
                                       TotMass;
                           atom1->z -= atom2->mass * MoveVector.z /
                                       TotMass;

                           atom2->x += atom1->mass * MoveVector.x /
                                       TotMass;
                           atom2->y += atom1->mass * MoveVector.y /
                                       TotMass;
                           atom2->z += atom1->mass * MoveVector.z /
                                       TotMass;
                        }
                        else if(atom1->Mobile)
                        {
                           atom1->x -= MoveVector.x;
                           atom1->y -= MoveVector.y;
                           atom1->z -= MoveVector.z;
                        }
                        else /* atom2 must be mobile                    */
                        {
                           atom2->x += MoveVector.x;
                           atom2->y += MoveVector.y;
                           atom2->z += MoveVector.z;
                        }
                     }
                  }
               }
            }
         }
      }

      if(Finished)
         return(NIter);
   }

#if (MAXVDWITER > 1)
   sprintf(gError,"Failed to relax VDW contacts in %d iterations",
           MAXVDWITER);
   StoreError("VDWRelax()",gError);
#endif

   return(0);
}
