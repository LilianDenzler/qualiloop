#define DISTAL
/*************************************************************************

   Program:    ECalc
   File:       energy.c
   
   Version:    V1.5.2
   Date:       05.03.21
   Function:   Calculate the energy of the structure
   
   Copyright:  (c) UCL, Prof. Andrew C. R. Martin 1994-2021
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
   V0.1   01.09.94 Original
   V0.2   30.09.94 Added checks on number of parameters
   V1.0   30.09.94 First release version
   V1.1   12.10.94 Added code for residue pseudo-energy
   V1.2   10.11.94 Fixed bugs in reading confs. for second chain
   V1.3   28.11.94 Fixed bug in display of atom numbers
   V1.3a  05.01.95 Added #ifdef'd SHOW_HBONDS
   V1.3b  31.01.95 Initialisation of variables for gcc
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Fixed for new version of GetWord()
                   Ignore Nter proline-CD/HT1 or HT2 in grid
   V1.5.1 07.01.21 Removed unused variables
   V1.5.2 05.03.21 Variable initializations to deal with missing atoms

*************************************************************************/
/* Includes
*/
#define ENERGY_MAIN

#include <math.h>
#include <sys/types.h>
#include <time.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/angle.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

#include "ecalc.h"

/************************************************************************/
/* Defines and macros
*/
#define TAYLORCUT 0.1
#define GRIDCHUNK 10
#define ELECT_CONST ((REAL)332.0716)

/* There seems to be some confusion over the standard naming for this   */
#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC CLOCKS_PER_SECOND
#endif

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
#include "energy.p"

/************************************************************************/
/*>REAL EBond(TOPOLOGY *topol)
   ---------------------------
   Calculate the energy of the bonds in the topology

   01.09.94 Original    By: ACRM
   05.09.94 Added 0.0 distance check and removed /2.0
   07.09.94 Atoms stored as pointers rather than offsets
   13.09.94 Added check on valid coordinates
   19.09.94 Added check on gNumBondParams
*/
REAL EBond(TOPOLOGY *topol)
{
   REAL TotE = (REAL)0.0,
        dist,
        DeltaLen;
   BOND *p;
   ATOM *atom1, *atom2;
   
   if(topol->NBonds == 0 || gNumBondParams == 0)
      return(TotE);

   for(p=topol->bonds; p!=NULL; NEXT(p))
   {
      /* atom1 & atom2 are pointers to type atom set by using offsets
         (p->atomX) into the topol->atom array
      */
      atom1     = p->atom1;
      atom2     = p->atom2;

      /* Check for valid coords                                         */
      if(VALID(atom1) && VALID(atom2))
      {
         if((dist  = DIST(atom1,atom2)) != 0.0)
         {
            DeltaLen  = dist - p->OptLen;
            TotE     += (DeltaLen * DeltaLen * p->force);
         }
      }
      
   }
   return(TotE);
}


/************************************************************************/
/*>REAL EAngle(TOPOLOGY *topol)
   ----------------------------
   Calculate the energy of the angles in the topology

   01.09.94 Original    By: ACRM
   05.09.94 Removed /2.0
   07.09.94 Atoms stored as pointers rather than offsets
   13.09.94 Added check on valid coordinates
   19.09.94 Added check on gNumAngleParams
*/
REAL EAngle(TOPOLOGY *topol)
{
   REAL  TotE = (REAL)0.0,
         ang,
         DeltaAng;
   ANGLE *p;
   ATOM  *atom1, *atom2, *atom3;
   
   if(topol->NAngles == 0 || gNumAngleParams == 0)
      return(TotE);

   for(p=topol->angles; p!=NULL; NEXT(p))
   {
      /* atom1, atom2 & atom3 are pointers to type atom set by using 
         offsets (p->atomX) into the topol->atom array
      */
      atom1     = p->atom1;
      atom2     = p->atom2;
      atom3     = p->atom3;

      if(VALID(atom1) && VALID(atom2) && VALID(atom3))
      {
         ang       = angle(atom1->x,atom1->y,atom1->z,
                           atom2->x,atom2->y,atom2->z,
                           atom3->x,atom3->y,atom3->z);
         DeltaAng  = ang - p->OptAng;
         TotE     += (DeltaAng * DeltaAng * p->force);
      }
   }
   return(TotE);
}
   

/************************************************************************/
/*>REAL ETorsion(TOPOLOGY *topol)
   ------------------------------
   Calculates the torsion energy.

   02.09.94 Original    By: ACRM
   07.09.94 Atoms stored as pointers rather than offsets
   13.09.94 Added check on valid coordinates
   15.09.94 Added chain name to error message
   19.09.94 Added check on gNumTorsionParams
   28.11.94 Calls AtomCount() to find atom numbers. 
            Added conf num to error message.
*/
REAL ETorsion(TOPOLOGY *topol)
{
   REAL     TotE = (REAL)0.0,
            energy,
            tor,
            force,
            CosTor,
            CosTorSq;
   TORSION  *p;
   ATOM     *atom1, *atom2, *atom3, *atom4;

   /* Check we've got some torsions                                     */
   if(topol->NTorsions == 0 || gNumTorsionParams == 0)
      return(TotE);

   /* For each torsion angle                                            */
   for(p=topol->torsions; p!=NULL; NEXT(p))
   {
      REAL ang1, ang2;
      
      /* atom1, atom2, atom3 & atom4 are pointers to type atom set by 
         using offsets (p->atomX) into the topol->atom array
      */
      atom1 = p->atom1;
      atom2 = p->atom2;
      atom3 = p->atom3;
      atom4 = p->atom4;

      if(VALID(atom1) && VALID(atom2) && VALID(atom3) && VALID(atom4))
      {
         /* Check for linearity                                         */
         ang1 = angle(atom1->x,atom1->y,atom1->z,
                      atom2->x,atom2->y,atom2->z,
                      atom3->x,atom3->y,atom3->z);
         ang2 = angle(atom2->x,atom2->y,atom2->z,
                      atom3->x,atom3->y,atom3->z,
                      atom4->x,atom4->y,atom4->z);
         if(ABS(ang1) <= (double)0.1 || ABS(ang2) <= (double)0.1)
         {
            sprintf(gError,"Warning: Conf. %d, Torsion almost linear \
(Chain %c Atoms %d %d %d %d)\n",
                    gConfNum,
                    topol->ChainName,
                    AtomCount(topol->atoms, atom1),
                    AtomCount(topol->atoms, atom2),
                    AtomCount(topol->atoms, atom3),
                    AtomCount(topol->atoms, atom4));
            StoreError("ETorsion()",gError);
         }
         
         /* Calculate torsion and derived values                        */
         tor       = phi(atom1->x,atom1->y,atom1->z,
                         atom2->x,atom2->y,atom2->z,
                         atom3->x,atom3->y,atom3->z,
                         atom4->x,atom4->y,atom4->z);
         CosTor    = (REAL)cos((double)tor);
         CosTorSq  = CosTor*CosTor;
         
         switch(p->period)
         {
         case 1:
            energy = CosTor;
            break;
         case 2:
            energy = (REAL)2.0*CosTorSq - (REAL)1.0;
            break;
         case 3:
            energy = CosTor * ((REAL)4.0*CosTorSq - (REAL)3.0);
            break;
         case 4:
            energy = (REAL)1.0 + CosTorSq*(REAL)8.0*(CosTorSq-(REAL)1.0);
            break;
         case 6:
            energy = CosTorSq*(CosTorSq*(CosTorSq*(REAL)32.0 - 
                                         (REAL)48.0) + 
                               (REAL)18.0) - (REAL)1.0;
            break;
         default:
            /* Case 5 is illegal and will have been screened out by
               StoreTorsionParam()
               */
            continue;
         }
         
         force = p->force;
         if(ABS(p->OptTor) > (REAL)0.01)
            force = (-force);
         
         energy = p->force + force*energy;
         
         TotE += energy;
      }
   }

   return(TotE);
}


/************************************************************************/
/*>REAL EImproper(TOPOLOGY *topol)
   -------------------------------
   Calculates the improper torsion energy.

   02.09.94 Original    By: ACRM
   07.09.94 Atoms stored as pointers rather than offsets
   13.09.94 Added check on valid coordinates
   15.09.94 Added chain name to error message
   19.09.94 Added check on gNumImproperParams
   28.11.94 Calls AtomCount() to find atom numbers
            Added conf num to error message.
*/
REAL EImproper(TOPOLOGY *topol)
{
   REAL     TotE  = (REAL)0.0,
            pi    = PI,
            TwoPi = (REAL)2.0 * PI,
            energy,
            imp,
            force,
            DeltaImp,
            SinImp;
   IMPROPER *p;
   ATOM     *atom1, *atom2, *atom3, *atom4;

   if(topol->NImpropers == 0 || gNumImproperParams == 0)
      return(TotE);

   for(p=topol->impropers; p!=NULL; NEXT(p))
   {
      REAL ang1, ang2;
      
      /* atom1, atom2, atom3 & atom4 are pointers to type atom set by 
         using offsets (p->atomX) into the topol->atom array
      */
      atom1 = p->atom1;
      atom2 = p->atom2;
      atom3 = p->atom3;
      atom4 = p->atom4;

      if(VALID(atom1) && VALID(atom2) && VALID(atom3) && VALID(atom4))
      {
         /* Check for linearity                                         */
         ang1 = angle(atom1->x,atom1->y,atom1->z,
                      atom2->x,atom2->y,atom2->z,
                      atom3->x,atom3->y,atom3->z);
         ang2 = angle(atom2->x,atom2->y,atom2->z,
                      atom3->x,atom3->y,atom3->z,
                      atom4->x,atom4->y,atom4->z);
         if(ABS(ang1) <= (double)0.1 || ABS(ang2) <= (double)0.1)
         {
            sprintf(gError,"Warning: Conf. %d, Improper almost linear \
(Chain %c Atoms %d %d %d %d)\n",
                    gConfNum,
                    topol->ChainName,
                    AtomCount(topol->atoms, atom1),
                    AtomCount(topol->atoms, atom2),
                    AtomCount(topol->atoms, atom3),
                    AtomCount(topol->atoms, atom4));
            StoreError("EImproper()",gError);
         }
         
         imp       = phi(atom1->x,atom1->y,atom1->z,
                         atom2->x,atom2->y,atom2->z,
                         atom3->x,atom3->y,atom3->z,
                         atom4->x,atom4->y,atom4->z);
         SinImp    = (REAL)sin((double)imp);
         
         DeltaImp = imp - p->OptTor;
         
         if(DeltaImp >  pi) DeltaImp -= TwoPi;
         if(DeltaImp < -pi) DeltaImp += TwoPi;
         
         if(p->OptTor == 0.0)
         {
            /* TAYLORCUT can be modified to minimise error              */
            if(ABS(SinImp) > TAYLORCUT)
            {
               force   = DeltaImp * p->force;
               energy  = DeltaImp * force;
            }
            else
            {
               energy  = imp * imp * p->force;
            }
         }
         else
         {
            if(ABS(SinImp) < (REAL)0.001) SinImp = (REAL)0.001;
            if(ABS(DeltaImp) >= pi/2.0)
            {
               sprintf(gError,"Warning: Conf. %d, Bent improper torsion \
too far from mimimum (Chain %c Atoms %d %d %d %d)",
                       gConfNum,
                       topol->ChainName,
                       AtomCount(topol->atoms, atom1),
                       AtomCount(topol->atoms, atom2),
                       AtomCount(topol->atoms, atom3),
                       AtomCount(topol->atoms, atom4));
            }
            force   = DeltaImp * p->force;
            energy  = DeltaImp * force;
         }
         TotE += energy;
      }
   }

   return(TotE);
}
   

/************************************************************************/
/*>GRID *BuildGrid(MOLECULE *mol, REAL GridCut)
   --------------------------------------------
   Allocate memory and build the non-bonded contacts grid

   02.09.94 Original    By: ACRM
   07.09.94 Ignores disulphide topology
   08.09.94 Modified for topol->atoms being an array of pointers
   13.09.94 Added check on valid coordinates
   15.09.94 Checks that neighbours have been allocated before freeing
   07.02.03 Added code to check for Nter Pro - if so don't include
            HT1/CD or HT2/CD in the grid
   05.03.21 Added initializations
*/
GRID *BuildGrid(MOLECULE *mol, REAL GridCut)
{
   int  AtomNum,
        Chain1,
        Chain2,
        i, j;
   REAL GridCutSq = GridCut * GridCut;
   BOOL NTerPro;
   
   /* If there is already a grid, free it                               */
   if(gGrid != NULL)
   {
      for(AtomNum=0; AtomNum<mol->NAtoms; AtomNum++)
      {
         if(gGrid[AtomNum].NeighbMax)
            free(gGrid[AtomNum].neighbours);
      }
      free(gGrid);
   }

   /* Allocate memory for the gGrid array                               */
   if((gGrid = (GRID *)malloc(mol->NAtoms * sizeof(GRID))) == NULL)
      return(NULL);
   for(i=0; i<mol->NAtoms; i++)
   {
      gGrid[i].NNeighb    = 0;
      gGrid[i].NeighbMax  = 0;
      /* 05.03.21 Added these intializations                            */
      gGrid[i].atom       = NULL;
      gGrid[i].antecedant = NULL;
      gGrid[i].neighbours = NULL;
   }
   
   AtomNum = 0;

   /* For each chain                                                    */
   for(Chain1=0; Chain1 < mol->NChains; Chain1++)
   {
      if(!(mol->topol[Chain1])->Disulphide)
      {
         /* For each atom in the chain                                  */
         for(i=0; i < (mol->topol[Chain1])->NAtoms; i++)
         {
            if(VALID(mol->topol[Chain1]->atoms[i]))
            {
               /* Put this atom in the grid                             */
               gGrid[AtomNum].atom = mol->topol[Chain1]->atoms[i];

               /* Check whether we have an Nterminal proline            */
               if((!strncmp(mol->topol[Chain1]->atoms[i]->atom,
                            "HT1 ", 4) ||
                   !strncmp(mol->topol[Chain1]->atoms[i]->atom,
                            "HT2 ", 4)) &&
                  mol->topol[Chain1]->sequence[1] == 'P')
               {
                  NTerPro = TRUE;
               }
               else
               {
                  NTerPro = FALSE;
               }

               
               /* For each remaining atom in the chain                  */
               for(j=(i+1); j < (mol->topol[Chain1])->NAtoms; j++)
               {
                  /* If we have NTer proline, ignore interactions
                     between Pro-CD and the NTer hydrogens
                  */
                  if(NTerPro && 
                     !strncmp(mol->topol[Chain1]->atoms[j]->atom,
                              "CD  ", 4))
                  {
                     continue;
                  }

                  if(!AddToGrid(gGrid, AtomNum, 
                                mol->topol[Chain1]->atoms[j], GridCutSq))
                  {
                     StoreError("BuildGrid()","No memory to add to grid");
                     return(FALSE);
                  }
               }
               
               /* For each remaining chain                              */
               for(Chain2=(Chain1+1); Chain2 < mol->NChains; Chain2++)
               {
                  if(!(mol->topol[Chain2])->Disulphide)
                  {
                     /* For each atom in this other chain               */
                     for(j=0; j < (mol->topol[Chain2])->NAtoms; j++)
                     {
                        if(!AddToGrid(gGrid,AtomNum,
                                      mol->topol[Chain2]->atoms[j], 
                                      GridCutSq))
                        {
                           StoreError("BuildGrid()",
                                      "No memory to add to grid");
                           return(FALSE);
                        }
                     }
                  }
               }
               
               AtomNum++;
            }
         }
      }
   }
         
   return(gGrid);
}


/************************************************************************/
/*>BOOL AddToGrid(GRID *grid, int AtomNum, ATOM *atom, REAL GridCutSq)
   -------------------------------------------------------------------
   Checks to see whether atom is within range of the base atom at this
   grid position. If so adds it to the grid array, allocating memory
   if necessary

   02.09.94 Original    By: ACRM
*/
BOOL AddToGrid(GRID *grid, int AtomNum, ATOM *atom, REAL GridCutSq)
{
   int i;
   
   /* Check the exclusions for these atom                               */
   for(i=0; i<(grid[AtomNum].atom)->nexcl; i++)
   {
      if((grid[AtomNum].atom)->excl[i] == atom)
         return(TRUE);
   }
   for(i=0; i<atom->nexcl; i++)
   {
      if(grid[AtomNum].atom == atom->excl[i])
         return(TRUE);
   }
   
   /* If atoms are within grid range                                    */
   if(DISTSQ(grid[AtomNum].atom, atom) <= GridCutSq)
   {
      /* Allocate more memory if necessary                              */
      if(++(grid[AtomNum].NNeighb) > grid[AtomNum].NeighbMax)
      {
         if(grid[AtomNum].NeighbMax == 0)
         {
            grid[AtomNum].NeighbMax = GRIDCHUNK;
            if((grid[AtomNum].neighbours = (ATOM **)
                malloc(grid[AtomNum].NeighbMax * sizeof(ATOM *)))==NULL)
            {
               StoreError("AddToGrid()","No memory for grid");
               return(FALSE);
            }
         }
         else
         {
            grid[AtomNum].NeighbMax += GRIDCHUNK;
            if((grid[AtomNum].neighbours = (ATOM **)
                realloc(grid[AtomNum].neighbours,
                        grid[AtomNum].NeighbMax * sizeof(ATOM *)))==NULL)
            {
               StoreError("AddToGrid()","No memory for grid");
               return(FALSE);
            }
         }
      }

      /* Store the atom pointer in the grid array                       */
      grid[AtomNum].neighbours[grid[AtomNum].NNeighb - 1] = atom;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void TrimGrid(int NGrid, GRID *grid, ZONE *zones, ZONE *ignores)
   ----------------------------------------------------------------
   Trims the non-bonded contacts grid to include only entries where
   either the base atom or one of the neighbours is in the specified
   zones. For entries included because a neighbour is in range, all
   neighbours not in range are removed.

   16.09.94 Original    By: ACRM
   21.09.94 Also removes ignored atoms
*/
void TrimGrid(int NGrid, GRID *grid, ZONE *zones, ZONE *ignores)
{
   int  i,
        j;
   BOOL OneInZone;

   /* If no zones or ignores have been defined, simply return           */
   if(zones == NULL && ignores == NULL)
      return;

   /* Run through each base atom in the grid                            */
   for(i=0; i<NGrid; i++)
   {
      /* First handle zones                                             */
      if(zones != NULL)
      {
         /* If the base atom is not in any of the zones                 */
         if(!InZones(grid[i].atom, zones))
         {
            OneInZone = FALSE;
            
            /* Run through the neighbours and see if any of these is in a
               zone, removing any which are not
            */
            for(j=0; j<grid[i].NNeighb; j++)
            {
               if(InZones(grid[i].neighbours[j], zones))
                  OneInZone = TRUE;
               else
                  grid[i].neighbours[j] = NULL;
            }
            
            /* If none of the neighbours was in zone, set this base atom 
               to NULL
            */
            if(!OneInZone)
               grid[i].atom = NULL;
         }
      }

      /* Now handle ignores                                             */
      if(ignores != NULL)
      {
         /* If the base atom is found in any of the zones, remove entry
            from grid
         */
         if(InZones(grid[i].atom, ignores))
         {
            grid[i].atom = NULL;
         }
         else
         {
            BOOL OneOK = FALSE;
            
            /* Run through the neighbours and see if any of these is to
               be ignored. If so, remove.
            */
            for(j=0; j<grid[i].NNeighb; j++)
            {
               if(InZones(grid[i].neighbours[j], ignores))
                  grid[i].neighbours[j] = NULL;
               else
                  OneOK = TRUE;
            }
            
            /* If all neighbours were in ignore zones, set this base 
               atom to NULL
            */
            if(!OneOK)
               grid[i].atom = NULL;
         }
      }
   }
}


/************************************************************************/
/*>REAL ENonBond(MOLECULE *mol, EPARAMS *eparams, 
                 REAL *VdwAE, REAL *VdwRE, REAL *ElectE)
   -----------------------------------------------------
   Calculates the non-bonded energy (electrostatic and Lennard-Jones)

   02.09.94 Original    By: ACRM
   07.09.94 Added check on zero distance
   13.09.94 Outputs individual energy components
   16.09.94 Checks that the atoms are valid
   19.09.94 Added check on gNumNonBondParams
*/
REAL ENonBond(MOLECULE *mol, EPARAMS *eparams, 
              REAL *VdwAE, REAL *VdwRE, REAL *ElectE)
{
   int  i,
        j;
   ATOM *atom1,
        *atom2;
   REAL CutNBSq,
        DistSq,
        InvDist2,
        InvDist4,
        InvDist6,
        InvDist12,
        energy,
        ParamR6,
        ParamR12;

   if(gNumNonBondParams == 0) return((REAL)0.0);

   *ElectE = (REAL)0.0;
   *VdwAE  = (REAL)0.0;
   *VdwRE  = (REAL)0.0;

   if(gGrid == NULL)
   {
      StoreError("ENonBond()","No non-bond grid has been defined");
      return((REAL)0.0);
   }

   CutNBSq = eparams->NonBondCut * eparams->NonBondCut;
   
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
               if((DistSq=DISTSQ(atom1, atom2)) < CutNBSq)
               {
                  if(DistSq != (REAL)0.0)
                  {
                     InvDist2 = (REAL)1.0 / DistSq;
                     
                     /* First the electrostatic energy                  */
                     if(eparams->ConstDielectric)
                     {
                        energy = atom1->charge * atom2->charge *
                           ELECT_CONST * (REAL)sqrt((double)InvDist2) / 
                              eparams->eta;
                     }
                     else
                     {
                        /* Distance dependent dielectric (eta = r)      */
                        energy = atom1->charge * atom2->charge * 
                           ELECT_CONST * InvDist2; 
                     }
                     *ElectE += energy;
                     
                     /* Now the Lennard-Jones energy                    */
                     InvDist4  = InvDist2 * InvDist2;
                     InvDist6  = InvDist2 * InvDist4;
                     InvDist12 = InvDist6 * InvDist6;
                     
                     /* Look up the parameters for this atom pair       */
                     ParamR12=gNonBondParams[atom1->key][atom2->key].r12;
                     ParamR6 =gNonBondParams[atom1->key][atom2->key].r6;

#ifdef DEBUG
                     {
                        REAL r,a;
                        r =  (ParamR12 * InvDist12);
                        a =  -1.0 * (ParamR6  * InvDist6);
                        fprintf(stderr, "Atoms: %s %s Attr %.6f \
Rep %.6f\n",
                                atom1->atom, atom2->atom, a, r);
                     }
#endif
                     
                     *VdwRE += (ParamR12 * InvDist12);
                     *VdwAE -= (ParamR6  * InvDist6);
                  }
               }
            }
         }
      }
   }

   return((*ElectE) + (*VdwRE) + (*VdwAE));
}


/************************************************************************/
/*>REAL EHBond(MOLECULE *mol, EPARAMS *eparams)
   --------------------------------------------
   Calculates the hydrogen bond energy

   05.09.94 Original    By: ACRM
   16.09.94 Checks atoms are valid
   19.09.94 Added check on gNumHBondParams
   05.01.95 Added #ifdef'd SHOW_HBONDS
   31.01.95 Initialised Rul3 and Rua3 for gcc
*/
REAL EHBond(MOLECULE *mol, EPARAMS *eparams)
{
   REAL CutOnHBSq,
        CutOffHBSq,
        CutOnHBAngSq,
        CutOffHBAngSq,
        CosAng,
        CosAngSq,
        ETot,
        EAng,
        DistSq,
        InvDistSq,
        InvDist10,
        energy,
        ParamR10,
        ParamR12,
        Rul3 = (REAL)0.0,
        /* Rul12,                                                       */
        Rua3 = (REAL)0.0;
        /* Rua12;                                                       */

   ATOM *AtomD,
        *AtomH,
        *AtomA;
   int  i,
        j;

   if(gNumHBondParams == 0) return((REAL)0.0);

   ETot     = (REAL)0.0;

   if(gHBGrid == NULL)
   {
      StoreError("EHBond()","No HBond grid has been defined");
      return(ETot);
   }

   CutOnHBSq  = eparams->CutOnHB  * eparams->CutOnHB;
   CutOffHBSq = eparams->CutOffHB * eparams->CutOffHB;

   if(CutOffHBSq != CutOnHBSq)
   {
      Rul3  = (REAL)1.0/((CutOffHBSq - CutOnHBSq) * 
                         (CutOffHBSq - CutOnHBSq) * 
                         (CutOffHBSq - CutOnHBSq));
      /* Rul12 = (REAL)12.0 * Rul3;                                     */
   }
   
   CutOnHBAngSq   = cos(eparams->CutOnHBAng);
   CutOnHBAngSq  *= CutOnHBAngSq;
   CutOffHBAngSq  = cos(eparams->CutOffHBAng);
   CutOffHBAngSq *= CutOffHBAngSq;

   if(CutOffHBAngSq != CutOnHBAngSq)
   {
      Rua3  = (REAL)1.0/((CutOffHBAngSq - CutOnHBAngSq) * 
                         (CutOffHBAngSq - CutOnHBAngSq) * 
                         (CutOffHBAngSq - CutOnHBAngSq));
      /* Rua12 = (REAL)12.0 * Rua3;                                     */
   }

   for(i=0; i<mol->NDonors; i++)
   {
      AtomH = gHBGrid[i].atom;
      AtomD = gHBGrid[i].antecedant;

      if(VALID(AtomH) && VALID(AtomD))
      {
         for(j=0; j<gHBGrid[i].NNeighb; j++)
         {
            AtomA = gHBGrid[i].neighbours[j];
            
            if(VALID(AtomA))
            {
               DistSq = DISTSQ(AtomD, AtomA);
               
               if(DistSq != (REAL)0.0)
               {
                  /* If we're within the HBond cutoff, calculate energy */
                  if(DistSq < CutOffHBSq)
                  {
                     InvDistSq = (REAL)1.0 / DistSq;
                     InvDist10 = InvDistSq * InvDistSq * InvDistSq * 
                                 InvDistSq * InvDistSq;
                     
                     /* Look up the parameters for this atom pair       */
                     ParamR12 = gHBondParams[AtomD->key][AtomA->key].r12;
                     ParamR10 = gHBondParams[AtomD->key][AtomA->key].r10;
                     
                     energy = (ParamR12 * InvDistSq * InvDist10) - 
                              (ParamR10 * InvDist10);
                     
                     /* If we're above the start of the smoothing range, 
                        calculate the smoothing factor.
                     */
                     if(DistSq > CutOnHBSq)
                     {
                        REAL DistFromOn,
                        DistFromOff,
                        Smoothing;
                        
                        DistFromOn  = CutOnHBSq  - DistSq;
                        DistFromOff = CutOffHBSq - DistSq;
                        
                        Smoothing = DistFromOff * DistFromOff * Rul3 *
                           (DistFromOff - (REAL)3.0 * DistFromOn);
                        
                        energy *= Smoothing;
                     }
                     
                     if(AtomH == NULL)
                     {
                        /* No explicit hydrogen, just add this onto the 
                           total energy
                        */
                        ETot += energy;
                     }
                     else
                     {
                        /* If we've got an explicit hydrogen, calculate 
                           the angle contribution
                        */
                        CosAng = (REAL)cos((double)
                                 angle(AtomD->x, AtomD->y, AtomD->z,
                                       AtomH->x, AtomH->y, AtomH->z,
                                       AtomA->x, AtomA->y, AtomA->z));
                        if(CosAng <= (REAL)(-0.99999))
                           CosAng = (REAL)(-0.99999);
                        
                        if(CosAng <= 0.0)
                        {
                           CosAngSq = CosAng * CosAng;
                           if(CosAngSq > CutOffHBAngSq)
                           {
                              EAng = CosAngSq * CosAngSq;
                              
                              if(CosAngSq < CutOnHBAngSq)
                              {
                                 REAL AngFromOn,
                                 AngFromOff,
                                 Smoothing;
                                 
                                 AngFromOn  = CutOnHBAngSq  - CosAngSq;
                                 AngFromOff = CutOffHBAngSq - CosAngSq;
                                 Smoothing  = AngFromOff * AngFromOff * 
                                    Rua3 *
                                    (AngFromOff - (REAL)3.0 * AngFromOn);
                                 
                                 EAng *= Smoothing;
                              }
                              
                              /* Add onto the total energy              */
                              ETot += (EAng * energy);
                           }
                        }
                     }
#ifdef SHOW_HBONDS
                     printf("Donor %4s Acceptor %4s DistDA %6.3f AngDHA \
%7.2f\n", AtomD->type, AtomA->type, 
          sqrt(DistSq), (180.0 * acos(CosAng) / PI)); 
#endif
                  }
               }
            }
         }
      }
   }

   return(ETot);
}


/************************************************************************/
/*>GRID *BuildHBGrid(MOLECULE *mol, REAL GridCut)
   ----------------------------------------------
   Allocate memory and build the H-bond contacts grid

   05.09.94 Original    By: ACRM
   13.09.94 Added check on valid coordinates
   15.09.94 Checks that neighbours have been allocated before freeing
   05.03.21 Added initializations for missing atoms
*/
GRID *BuildHBGrid(MOLECULE *mol, REAL GridCut)
{
   int      AtomNum,
            Chain1,
            Chain2,
            i;
   REAL     GridCutSq = GridCut * GridCut;
   DONOR    *d;
   ACCEPTOR *a;
   
   /* If there is already a grid, free it                               */
   if(gHBGrid != NULL)
   {
      for(AtomNum=0; AtomNum<mol->NDonors; AtomNum++)
      {
         if(gHBGrid[AtomNum].NeighbMax) 
            free(gHBGrid[AtomNum].neighbours);
      }
      free(gHBGrid);
   }

   /* Allocate memory for the gHBGrid array                             */
   if((gHBGrid = (GRID *)malloc(mol->NDonors * sizeof(GRID))) == NULL)
      return(NULL);
   for(i=0; i<mol->NDonors; i++)
   {
      gHBGrid[i].NNeighb    = 0;
      gHBGrid[i].NeighbMax  = 0;
      /* 05.03.21 Added these intializations                            */
      gHBGrid[i].atom       = NULL;
      gHBGrid[i].antecedant = NULL;
      gHBGrid[i].neighbours = NULL;
   }
   
   AtomNum = 0;

   /* For each chain                                                    */
   for(Chain1=0; Chain1 < mol->NChains; Chain1++)
   {
      if(!(mol->topol[Chain1])->Disulphide)
      {
         /* For each donor in the chain                                 */
         for(d=mol->topol[Chain1]->donors; d!=NULL; NEXT(d))
         {
            if(VALID(d->atom1) && VALID(d->atom2))
            {
               /* Put this donor atom in the grid                       */
               gHBGrid[AtomNum].atom       = d->atom1;
               
               /* Put the antecedant in the grid                        */
               gHBGrid[AtomNum].antecedant = d->atom2;
               
               /* For each chain                                        */
               for(Chain2=0; Chain2 < mol->NChains; Chain2++)
               {
                  if(!(mol->topol[Chain2])->Disulphide)
                  {
                     /* For each acceptor in the chain                  */
                     for(a=mol->topol[Chain2]->acceptors; 
                         a!=NULL; 
                         NEXT(a))
                     {
                        if(!AddToGrid(gHBGrid, AtomNum, 
                                      a->atom, 
                                      GridCutSq))
                        {
                           StoreError("BuildHBGrid()",
                                      "No memory to add to grid");
                           return(FALSE);
                        }
                     }
                  }
               }
               
               AtomNum++;
            }
         }
      }
   }
         
   return(gHBGrid);
}

/************************************************************************/
/*>BOOL DoEnergyCalculations(FILE *out, MOLECULE *mol, EPARAMS EParams, 
                             FLAGS Flags, FILE *conffp)
   --------------------------------------------------------------------
   Main routine to read conformations, calculate and print energies

   13.09.94 Original    By: ACRM
   14.09.94 Removed unused variable
   15.09.94 Added re-gridding & timings
   16.09.94 Added calls to TrimGrid()
            Build grid on Flags.elect
   21.09.94 Updated calls to TrimGrid() for ignores
   28.11.94 Uses gConfNum rather than ConfNum
*/
BOOL DoEnergyCalculations(FILE *out, MOLECULE *mol, EPARAMS EParams, 
                          FLAGS Flags, FILE *conffp)
{
   REAL    energy;
   clock_t time1,
           NBGridTime = 0,
           HBGridTime = 0;

   gConfNum = 0;

   if(Flags.vdwa || Flags.vdwr || Flags.elect)
   {
      /* If there is no close contacts grid, build one                  */
      if(gGrid == NULL)
      {
         time1 = clock();
      
         if((gGrid = BuildGrid(mol, EParams.GridCut))==NULL)
         {
            StoreError("DoEnergyCalculations()","Unable to build grid");
            return(FALSE);
         }
         TrimGrid(mol->NAtoms, gGrid, gZone, gIgnore);
         
         if(Flags.Timings)
            NBGridTime = clock() - time1;
      }
   }
   
   if(Flags.hbonds)
   {
      /* If there is no HBond contacts grid, build one                  */
      if(gHBGrid == NULL)
      {
         time1 = clock();

         if((gHBGrid = BuildHBGrid(mol, EParams.GridCut))==NULL)
         {
            StoreError("DoEnergyCalculations()","Unable to build HBond \
grid");
            return(FALSE);
         }
         TrimGrid(mol->NDonors, gHBGrid, gZone, gIgnore);
         
         if(Flags.Timings)
            HBGridTime = clock() - time1;
      }
   }

   if(Flags.Timings)
      fprintf(out, "\nTIMES: Non-bond grid %5.2f, H-bond grid %5.2f\n",
              (REAL)NBGridTime/(REAL)CLOCKS_PER_SEC,
              (REAL)HBGridTime/(REAL)CLOCKS_PER_SEC);
   
   /* Calculate the base structure energy                               */
   energy = ShowEnergy(out, mol, &EParams, &Flags);

   fprintf(out,"Base structure energy = %f%s\n",energy,
           gAtomError?" (MISSING ATOMS)":"");

   /* For each conformation                                             */
   while(GetConf(conffp, mol, TRUE))
   {
      gConfNum++;
      
      if(Flags.Regrid && !(gConfNum % Flags.Regrid))
      {
         if(Flags.vdwa || Flags.vdwr || Flags.elect)
         {
            time1 = clock();
      
            if((gGrid = BuildGrid(mol, EParams.GridCut))==NULL)
            {
               StoreError("DoEnergyCalculations()",
                          "Unable to rebuild grid");
               return(FALSE);
            }
            TrimGrid(mol->NAtoms, gGrid, gZone, gIgnore);

            if(Flags.Timings)
               NBGridTime = clock() - time1;
         }
         if(Flags.hbonds)
         {
            time1 = clock();
      
            if((gHBGrid = BuildHBGrid(mol, EParams.GridCut))==NULL)
            {
               StoreError("DoEnergyCalculations()","Unable to build \
HBond grid");
               return(FALSE);
            }
            TrimGrid(mol->NDonors, gHBGrid, gZone, gIgnore);

            if(Flags.Timings)
               HBGridTime = clock() - time1;
         }

         if(Flags.Timings)
            fprintf(out, "\nTIMES: Non-bond grid %5.2f, \
H-bond grid %5.2f\n",
                    (REAL)NBGridTime/(REAL)CLOCKS_PER_SEC,
                    (REAL)HBGridTime/(REAL)CLOCKS_PER_SEC);
      }

      energy = ShowEnergy(out, mol, &EParams, &Flags);
      fprintf(out,"Conformation %d energy = %f\n", gConfNum, energy);

      CacheConf(gConfNum, energy, Flags.NCache);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>REAL ShowEnergy(FILE *out, MOLECULE *mol, EPARAMS *eparams, 
                   FLAGS *flags)
   -----------------------------------------------------------
   Main energy display routine

   13.09.94 Original    By: ACRM
   15.09.94 Added time code
   11.10.94 Added call to EResidue() and other residue code
   19.05.95 Added call to RelaxStructure()
*/
REAL ShowEnergy(FILE *out, MOLECULE *mol, EPARAMS *eparams, FLAGS *flags)
{
   REAL    BondE     = (REAL)0.0,
           AngleE    = (REAL)0.0,
           TorsionE  = (REAL)0.0,
           ImproperE = (REAL)0.0,
           VdwAE     = (REAL)0.0,
           VdwRE     = (REAL)0.0,
           HBondE    = (REAL)0.0,
           ElectE    = (REAL)0.0,
           /* NonBondE  = (REAL)0.0,                                    */
           ResidueE  = (REAL)0.0,
           ETot      = (REAL)0.0;
   int     i;
   clock_t TimeBonds     = 0,
           TimeAngles    = 0,
           TimeTorsions  = 0,
           TimeImpropers = 0,
           TimeNonBonds  = 0,
           TimeHBonds    = 0,
           TimeResidue   = 0,
           time1;
   ZONE    *zone;


   /* Relax the structure if required                                   */
   if(eparams->DoRelax)
      RelaxStructure(mol, eparams->VDWTol, eparams->ShakeTol);

   /* For each chain, calculate each energy component as required       */
   for(i=0; i<mol->NChains; i++)
   {
      if(flags->bonds)
      {
         time1 = clock();
         BondE     += EBond(mol->topol[i]);

         if(flags->Timings)
            TimeBonds += (clock()-time1);
      }
      
      if(flags->angles)
      {
         time1 = clock();
         AngleE    += EAngle(mol->topol[i]);

         if(flags->Timings)
            TimeAngles += (clock()-time1);
      }
      
      if(flags->torsions)  
      {
         time1 = clock();
         TorsionE  += ETorsion(mol->topol[i]);

         if(flags->Timings)
            TimeTorsions += (clock()-time1);
      }
      
      if(flags->impropers)
      {
         time1 = clock();
         ImproperE += EImproper(mol->topol[i]);

         if(flags->Timings)
            TimeImpropers += (clock()-time1);
      }
      
   }

   /* Do the non-bond energy                                            */
   if(flags->vdwa || flags->vdwr || flags->elect)
   {
      time1 = clock();
      /* NonBondE =                                                     */
      ENonBond(mol, eparams, &VdwAE, &VdwRE, &ElectE);

      if(flags->Timings)
         TimeNonBonds = (clock()-time1);

      /* If any of the non-bond flags were false, set the appropriate 
         values to zero
      */
      if(!flags->vdwa)  VdwAE  = (REAL)0.0;
      if(!flags->vdwr)  VdwRE  = (REAL)0.0;
      if(!flags->elect) ElectE = (REAL)0.0;
   }
   
   /* Do the H-bond energy                                              */
   if(flags->hbonds) 
   {
      time1 = clock();
      HBondE = EHBond(mol, eparams);

      if(flags->Timings)
         TimeHBonds += (clock()-time1);
   }

   /* Do the residue pseudo-energy                                      */
   if(flags->residue)
   {
      for(zone = gZone; zone!=NULL; NEXT(zone))
      {
         time1 = clock();
         ResidueE += EResidue(zone);

         if(flags->Timings)
            TimeResidue += (clock()-time1);
      }
   }
   
   /* Calculate total energy                                            */
   ETot = flags->BondScale     * BondE     + 
          flags->AngleScale    * AngleE    + 
          flags->TorsionScale  * TorsionE  + 
          flags->ImproperScale * ImproperE + 
          flags->VDWAScale     * VdwAE     + 
          flags->VDWRScale     * VdwRE     + 
          flags->HBondScale    * HBondE    + 
          flags->ResidueScale  * ResidueE  + 
          flags->ElectScale    * ElectE;

   if(out != NULL)   /* Display energy components                       */
   {
      if(flags->ShowBonds)     fprintf(out,"Bonds: %f ",          BondE);
      if(flags->ShowAngles)    fprintf(out,"Angles: %f ",         AngleE);
      if(flags->ShowTorsions)  fprintf(out,"Torsions: %f ",  TorsionE);
      if(flags->ShowImpropers) fprintf(out,"Impropers: %f ", ImproperE);
      if(flags->ShowVDWA)      fprintf(out,"VDW-Attraction: %f ", VdwAE);
      if(flags->ShowVDWR)      fprintf(out,"VDW-Repulsion: %f ",  VdwRE);
      if(flags->ShowElect)     fprintf(out,"Electrostatics: %f ", ElectE);
      if(flags->ShowHBonds)    fprintf(out,"HBonds: %f ",         HBondE);
      if(flags->ShowResidue)   fprintf(out,"Residue: %f ",   ResidueE);

      if(flags->ShowBonds     || 
         flags->ShowAngles    || 
         flags->ShowTorsions  || 
         flags->ShowImpropers || 
         flags->ShowVDWA      || 
         flags->ShowVDWR      || 
         flags->ShowElect     || 
         flags->ShowResidue   || 
         flags->ShowHBonds)
      {
         fprintf(out,"\n");
      }

      if(flags->Timings)
      {
         fprintf(out,"TIMES: Bonds %4.2f, Angles %4.2f, Torsions %4.2f, \
Impropers %4.2f, Non-bonds %5.2f, H-bonds %5.2f, Residue %4.2f\n",
                 (REAL)TimeBonds     / (REAL)CLOCKS_PER_SEC,
                 (REAL)TimeAngles    / (REAL)CLOCKS_PER_SEC,
                 (REAL)TimeTorsions  / (REAL)CLOCKS_PER_SEC,
                 (REAL)TimeImpropers / (REAL)CLOCKS_PER_SEC,
                 (REAL)TimeNonBonds  / (REAL)CLOCKS_PER_SEC,
                 (REAL)TimeHBonds    / (REAL)CLOCKS_PER_SEC,
                 (REAL)TimeResidue   / (REAL)CLOCKS_PER_SEC);
      }
   }

   return(ETot);
}


/************************************************************************/
/*>BOOL GetConf(FILE *conffp, MOLECULE *mol, BOOL SkipReference)
   -------------------------------------------------------------
   Reads a conformation from a conf file into the molecule structure.
   On the first call, the header will be read from the file and 
   various checks will be made. If the SkipReference flag is the
   reference set of coordinates will be read on the first call and
   ignored.

   15.09.94 Original    By: ACRM
*/
BOOL GetConf(FILE *conffp, MOLECULE *mol, BOOL SkipReference)
{
   static ATOM **ConsAtoms = NULL;
   static int  ncons;
   int         i;
   char        buffer[MAXBUFF];
   
   if(conffp == NULL) return(FALSE);

   /* On the first call, read the header out of the conformation file   */
   if(ConsAtoms == NULL)
   {
      if((ConsAtoms = ReadConfHeader(conffp, mol, SkipReference, &ncons))
         ==NULL)
      {
         StoreError("GetConf()",
                    "Unable to read header from conformation file");
         return(FALSE);
      }
   }

   /* Read the next ncons records and store the coordinates             */
   for(i=0; i<ncons; i++)
   {
      if(!fgets(buffer,MAXBUFF,conffp))
         return(FALSE);
      sscanf(buffer, "%lf %lf %lf", &(ConsAtoms[i]->x),
                                    &(ConsAtoms[i]->y),
                                    &(ConsAtoms[i]->z));
   }

   /* Read the energy record                                            */
   fgets(buffer,MAXBUFF,conffp);

   return(TRUE);
}

/************************************************************************/
/*>ATOM **ReadConfHeader(FILE *conffp, MOLECULE *mol, BOOL SkipReference,
                         int *ncons)
   ----------------------------------------------------------------------
   Read the header out of a CSearch (Charmm-free CONGEN) conformation 
   file. Builds an array of atom pointers (NULL terminated) for the
   constructed atoms which is returned by the routine. The number of
   constructed atoms is output.

   If the SkipReference flag is set, the first conformation will also be
   read and ignored.

   15.09.94 Original    By: ACRM
   10.11.94 Fixed bug when constructed region not in first chain
            Added a missing NULL return on failure
   06.02.03 Fixed for new version of GetWord()
*/
ATOM **ReadConfHeader(FILE *conffp, MOLECULE *mol, BOOL SkipReference,
                      int *ncons)
{
   TOPOLOGY *topol;
   int      MolNAtoms,
            MolNRes,
            i, 
            nres, natoms, 
            NRead, 
            chain, 
            /* atomnum,                                                 */
            FoundAt,
            *ConsAtNum;
   ATOM     **ConsAtoms;
   char     buffer[MAXBUFF],
            word[16],
            *ch;
   
   /* Get the required number of residues atoms out of the MOLECULE
      structure
   */
   MolNAtoms = mol->NAtoms;
   MolNRes   = 0;
   for(i=0; i<mol->NChains; i++)
   {
      if(!mol->topol[i]->Disulphide)
         MolNRes += mol->topol[i]->NRes;
   }

   /* Get the number of residues, number of atoms and number of
      constructed atoms out of the file
   */
   if(fscanf(conffp,"%d %d %d",&nres, &natoms, ncons) != 3)
   {
      StoreError("ReadConfHeader()",
                 "Conformation file has invalid header");
      return(NULL);
   }
   
   /* Check the atom and residue counts match                           */
   if(MolNRes != nres || MolNAtoms != natoms)
   {
      StoreError("ReadConfHeader()",
                 "Conformation file does not match this structure");
      sprintf(gError, "Number or residues in topology: %d, in \
conformation file %d",MolNRes, nres);
      StoreError("                ",gError);
      sprintf(gError, "Number or atoms in topology:    %d, in \
conformation file %d",MolNAtoms, natoms);
      StoreError("                ",gError);

      return(NULL);
   }
   
   /* Allocate memory for the constructed atom numbers array            */
   if((ConsAtNum = (int *)malloc((*ncons) * sizeof(int))) == NULL)
   {
      StoreError("ReadConfHeader()",
                 "No memory for constructed atom number array");
      return(NULL);
   }
   
   /* Read the constructed atom numbers list                            */
   NRead = 0;
   while(NRead < (*ncons))
   {
      fgets(buffer,MAXBUFF,conffp);
      TERMINATE(buffer);
      
      ch = buffer; 
      while(ch != NULL && *ch)
      {
         ch = GetWord(ch, word, 16);
         if(sscanf(word,"%d",&(ConsAtNum[NRead])) == 0)
         {
            StoreError("ReadConfHeader()",
                       "Corrupted constructed atoms list in conformation \
file");
            return(NULL);
         }
         
         if(++NRead >= (*ncons))
            break;
      }
   }

   /* Allocate memory for the array of atom pointers                    */
   if((ConsAtoms = (ATOM **)malloc((*ncons+1) * sizeof(ATOM *)))==NULL)
   {
      free(ConsAtNum);
      StoreError("ReadConfHeader()",
                 "No memory for constructed atom array");
      return(NULL);   /* Added 10.11.94                                 */
   }
   
   /* Convert the constructed atom numbers to a list of atom pointers   */
   for(i=0; i<(*ncons); i++)
   {
      chain   = 0;
      /* atomnum = 0;                                                   */
      FoundAt = ConsAtNum[i];

      for(chain = 0; chain < mol->NChains; chain++)
      {
         topol = mol->topol[chain];
         if(!topol->Disulphide)
         {
            /* See if the constructed atom is in this chain 
               10.11.94 Fixed from mol->NAtoms
            */
            if(ConsAtNum[i] <= topol->NAtoms)
            {
               /* Yes! Get the atom pointer                             */
               ConsAtoms[i] = topol->atoms[(ConsAtNum[i] - 1)];
               FoundAt      = 0;
               break;
            }
            else      /* Atom not in this chain                         */
            {
               /* Decr. atom number by number of atoms in this chain    
                  10.11.94 Corrected from ConsAtoms[i]
               */
               ConsAtNum[i] -= topol->NAtoms;
            }
         }
      }

      /* Check that it was found in one of the chains                   */
      if(FoundAt)
      {
         sprintf(gError,"Constructed atom %d not in any chain",FoundAt);
         StoreError("ReadConfHeader()", gError);

         free(ConsAtNum);
         free(ConsAtoms);
         
         return(NULL);
      }
   }
   
   /* Put a NULL at the end of the ConsAtom array to terminate it       */
   ConsAtoms[*ncons] = NULL;

   /* Free the memory for the constructed atom numbers list             */
   free(ConsAtNum);

   /* Skip over the first (reference) conformation if the SkipReference
      flag is set
   */
   if(SkipReference)
   {
      /* Read *ncons lines out of the file for the coordinates and one
         more for the energy
      */
      for(i=0; i<=(*ncons); i++)
         fgets(buffer,MAXBUFF,conffp);
   }

   return(ConsAtoms);
}  


/************************************************************************/
/*>REAL EResidue(ZONE *zone)
   -------------------------
   Calculates a statistical pseudo-energy component for the specified
   zone. PatchZones() must have been called first so that the zone
   structure is completed.

   11.10.94 Original    By: ACRM
   12.10.94 Changed pi from static since SG C-compiler doesn't seem to
            like the function (PI) in a static declaration
*/
REAL EResidue(ZONE *zone)
{
   VEC3F    base;
   ATOM     *FirstCA = NULL,
            *LastCA  = NULL,
            *distal,
            *atom;
   TOPOLOGY *topol;
   REAL     SCAng,
            Energy = (REAL)0.0;
   int      resnum,
            i;
   REAL     pi = PI;
   
   /* Check a zone has been specified                                   */
   if(zone==NULL)
      return((REAL)0.0);

   /* Get the topology                                                  */
   topol = zone->topol;

   /* Find the first and last CAs in the zone                           */
   for(i=zone->atom1; i<zone->atom2; i++)
   {
      atom = topol->atoms[i];
      if(!strncmp(atom->atom,"CA  ",4))
      {
         if(FirstCA == NULL)
            FirstCA = atom;
         LastCA = atom;
      }
   }
   
   /* Find the midpoint between the ends of the zone                    */
   base.x = FirstCA->x + ((LastCA->x - FirstCA->x) / (REAL)2.0);
   base.y = FirstCA->y + ((LastCA->y - FirstCA->y) / (REAL)2.0);
   base.z = FirstCA->z + ((LastCA->z - FirstCA->z) / (REAL)2.0);

   /* Run through the zone, finding the CA and the distal atom of each
      residue
   */
   for(i=zone->atom1; i<zone->atom2; i++)
   {
      atom = topol->atoms[i];

      /* See if we've hit a C-alpha                                     */
      if(!strncmp(atom->atom,"CA  ",4))
      {
         /* Find the residue number for this atom                       */
         if((resnum = GetResFromAtomNum(topol, i)) != (-1))
         {
#ifdef DISTAL
            /* Identify the distal atom                                 */
            if((distal = FindDistalAtom(topol, 
                                        topol->sequence[resnum], 
                                        i, 
                                        topol->ResStart[resnum+1]))!=NULL)
            {
               /* Calculate the angle between the base, the CA and the
                  distal atom
               */
               if(distal == atom)
               {
                  /* If this distal atom was the CA, set angle to 90    */
                  SCAng = pi / (REAL)2.0;
               }
               else
               {
                  SCAng = angle(base.x,    base.y,    base.z, 
                                atom->x,   atom->y,   atom->z,
                                distal->x, distal->y, distal->z);
               }

               /* Calculate the actual pseudo-energy for this residue
                  and angle
               */

               Energy += ScoreResidue(topol->sequence[resnum], SCAng);
            }
#else
            if(topol->sequence[resnum] == 'G')
            {
               SCAng = pi / (REAL)2.0;
            }
            else
            {
               VEC3F SCCofG;
            
               /* Find CofG of sidechain                                */
               FindSidechainCG(topol, topol->ResStart[resnum],
                               topol->ResStart[resnum+1], &SCCofG);
               SCAng = angle(base.x,    base.y,    base.z, 
                             atom->x,   atom->y,   atom->z,
                             SCCofG.x,  SCCofG.y,  SCCofG.z);
            }

            Energy += ScoreResidue(topol->sequence[resnum], SCAng);
#endif
         }
      }
   }

   return(Energy);
}


/************************************************************************/
/*>REAL ScoreResidue(char restype, REAL angle)
   -------------------------------------------
   Calculates an actual pseudo-energy value for a given residue type and
   an angle.

   11.10.94 Original    By: ACRM
*/
REAL ScoreResidue(char restype, REAL angle)
{
   RESPARAM *p;
   int      box;
   REAL     energy;
   
   for(p=gResParam; p!=NULL; NEXT(p))
   {
      /* If we've found the residue in the gResParam linked list        */
      if(p->restype == restype)
      {
         /* Find which box this angle falls into                        */
         box = (int)(p->nbox * angle / p->maxval);
         while(box >= p->nbox)
            box--;
         
         /* Return the score for this box. This is:

                    /        1             \
            scale * | ---------------  - 1 |
                    \ max(value, .01)      /

            i.e. if value is 1, this will return 0
                 if value <= .01, this will return 99*scale

         */

         energy = (((REAL)1.0 / MAX(p->values[box], 0.01)) - (REAL)1.0) * 
                p->scale;

         /* Return the score for this box. This is:

            scale * 100.0 * ( 1 - value )

            i.e. if value is 1, this will return 0
                 if value is 0, this will return 100*scale
         */
/*
//         energy += (REAL)(100.0 * p->scale * (1 - p->values[box]));
*/

         return(energy);
      }
   }

   return((REAL)0.0);
}


/************************************************************************/
/*>ATOM *FindDistalAtom(TOPOLOGY *topol, char restype, int start, 
                        int stop)
   --------------------------------------------------------------
   Returns a pointer to the distal atom of the sidechain.

   11.10.94 Original    By: ACRM
*/
ATOM *FindDistalAtom(TOPOLOGY *topol, char restype, int start, int stop)
{
   char     disttype[8];
   RESPARAM *p;
   ATOM     *atom;
   int      i;
   
   disttype[0] = '\0';

   /* Find the required distal atom type for this residue               */
   for(p=gResParam; p!=NULL; NEXT(p))
   {
      if(p->restype == restype)
      {
         strcpy(disttype, p->distal);
         break;
      }
   }
   if(!disttype[0])
      return(NULL);

   /* Run through the atoms till we find the required type              */
   for(i=start; i<stop; i++)
   {
      atom = topol->atoms[i];
      if(!strncmp(atom->atom, disttype, 4))
         return(atom);
   }
   
   return(NULL);
}

/************************************************************************/
/*>void FindSidechainCG(TOPOLOGY *topol, int start, int stop, VEC3F *CofG)
   -----------------------------------------------------------------------
   Finds the centre of geometry of a sidechain (or of a range of 
   sidechains). The search runs from atom `start' in the sepcified
   topology to less than atom `stop'. N, CA, C, O and hydrogen atoms
   are ignored in the calculation. The result is placed in CofG.

   12.10.94 Original    By: ACRM
*/
void FindSidechainCG(TOPOLOGY *topol, int start, int stop, VEC3F *CofG)
{
   ATOM *a;
   int  i,
        count = 0;   

   CofG->x = CofG->y = CofG->z = (REAL)0.0;
   
   for(i=start; i<stop; i++)
   {
      a = topol->atoms[i];
      
      if(strncmp(a->atom, "N   ",4) &&
         strncmp(a->atom, "CA  ",4) &&
         strncmp(a->atom, "C   ",4) &&
         strncmp(a->atom, "O   ",4) &&
         a->atom[0] != 'H')
      {
         CofG->x += a->x;
         CofG->y += a->y;
         CofG->z += a->z;
         count++;
      }
   }
   
   if(count)
   {
      CofG->x /= count;
      CofG->y /= count;
      CofG->z /= count;
   }
}


/************************************************************************/
/*>int AtomCount(ATOM **AtomArray, ATOM *atom)
   -------------------------------------------
   Counts through the atom array until we find the required pointer.
   Returns the atom number for the specified atom

   28.11.94 Original    By: ACRM
*/
int AtomCount(ATOM **AtomArray, ATOM *atom)
{
   int i;

   for(i=0; AtomArray[i] <= atom; i++) ;
   return(i);
}

