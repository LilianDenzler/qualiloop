/*************************************************************************

   Program:    ECalc
   File:       PrintParams.c
   
   Version:    V1.5.2
   Date:       05.03.21
   Function:   Print the contents of the parameters arrays
   
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
   V0.1   26.08.94 Original
   V0.2   30.09.94 Added checks on number of parameters
   V1.0   30.09.94 First release version
   V1.1            Skipped
   V1.2            Skipped
   V1.3   28.11.94 Bug fixes in energy.c
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Skipped
   V1.5.1 07.01.21 Skipped
   V1.5.2 05.03.21 Skipped

*************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "params.h"

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
/*>void PrintParams(FILE *fp)
   --------------------------
   Print the contents of the parameters arrays

   26.08.94 Original    By: ACRM
   05.09.94 Modified for 2D gNonBondParams array
   13.09.94 Fixed bug in printing HBond array
   30.09.94 Added checks on number of params
*/
void PrintParams(FILE *fp)
{
   int i,
       j;

   fprintf(fp,"*** Contents of the Parameter file ***\n=================\
=====================\n\n");
   fprintf(fp,"There were %d bond parameters\n",     gNumBondParams);
   fprintf(fp,"There were %d angle parameters\n",    gNumAngleParams);
   fprintf(fp,"There were %d torsion parameters\n",  gNumTorsionParams);
   fprintf(fp,"There were %d improper parameters\n", gNumImproperParams);
   fprintf(fp,"There were %d non-bond parameters\n", gNumNonBondParams);
   fprintf(fp,"There were %d H-bond parameters\n\n", gNumHBondParams);
   
   if(gNumBondParams)
   {
      fprintf(fp,"Bond parameters\n");
      fprintf(fp,"Atom 1    Atom2        force            B0\n");
      fprintf(fp,"------------------------------------------\n");
      for(i=0; i<gNumBondParams; i++)
         fprintf(fp,"%s      %s      %8.3f      %8.3f\n",
                 gBondParams[i].atom1,
                 gBondParams[i].atom2,
                 gBondParams[i].force,
                 gBondParams[i].OptLen);
      fprintf(fp,"\n");
   }
   
   if(gNumAngleParams)
   {
      fprintf(fp,"Angle parameters\n");
      fprintf(fp,"Atom 1    Atom 2    Atom 3       \
force            A0\n");
      fprintf(fp,"---------------------------------\
-------------------\n");
      for(i=0; i<gNumAngleParams; i++)
         fprintf(fp,"%s      %s      %s      %8.3f      %8.3f\n",
                 gAngleParams[i].atom1,
                 gAngleParams[i].atom2,
                 gAngleParams[i].atom3,
                 gAngleParams[i].force,
                 gAngleParams[i].OptAng);
      fprintf(fp,"\n");
   }
   
   if(gNumTorsionParams)
   {
      fprintf(fp,"Torsion parameters\n");
      fprintf(fp,"Atom 1    Atom 2    Atom 3    Atom 4       \
force        period            T0\n");
      fprintf(fp,"-------------------------------------------\
---------------------------------\n");
      for(i=0; i<gNumTorsionParams; i++)
         fprintf(fp,"%s      %s      %s      %s      %8.3f      \
%8.3f      %8.3f\n",
                 gTorsionParams[i].atom1,
                 gTorsionParams[i].atom2,
                 gTorsionParams[i].atom3,
                 gTorsionParams[i].atom4,
                 gTorsionParams[i].force,
                 gTorsionParams[i].period,
                 gTorsionParams[i].OptTor);
      fprintf(fp,"\n");
   }
   
   if(gNumImproperParams)
   {
      fprintf(fp,"Improper parameters\n");
      fprintf(fp,"Atom 1    Atom 2    Atom 3    Atom 4       \
force            I0\n");
      fprintf(fp,"-------------------------------------------\
-------------------\n");
      for(i=0; i<gNumImproperParams; i++)
         fprintf(fp,"%s      %s      %s      %s      %8.3f      %8.3f\n",
                 gImproperParams[i].atom1,
                 gImproperParams[i].atom2,
                 gImproperParams[i].atom3,
                 gImproperParams[i].atom4,
                 gImproperParams[i].force,
                 gImproperParams[i].OptTor);
      fprintf(fp,"\n");
   }
   
   if(gNumNonBondParams)
   {
      fprintf(fp,"Non-Bond parameters\n");
      fprintf(fp,"Atom       Polariz          NEff          \
VDWR          Mass\n");
      fprintf(fp,"------------------------------------------\
------------------\n");
      for(i=0; i<gNumAtomParams; i++)
         fprintf(fp,"%s      %8.3f      %8.3f      %8.3f      %8.3f\n",
                 gAtomParams[i].atom,
                 gAtomParams[i].pol,
                 gAtomParams[i].NEff,
                 gAtomParams[i].vdwr,
                 gAtomParams[i].mass);
      fprintf(fp,"\n");
      
      fprintf(fp,"Atom 1    Atom 2          R6           R12\n");
      fprintf(fp,"------------------------------------------\n");
      for(i=0; i<gNumAtomParams; i++)
      {
         for(j=i; j<gNumAtomParams; j++)
         {
            fprintf(fp,"%s      %s      %8.3f   %11.3f\n",
                    gNonBondParams[i][j].atom1,
                    gNonBondParams[i][j].atom2,
                    gNonBondParams[i][j].r6,
                    gNonBondParams[i][j].r12);
         }
      }
      fprintf(fp,"\n");
   }
   
   if(gNumHBondParams)
   {
      fprintf(fp,"H-Bond parameters\n");
      fprintf(fp,"Atom 1    Atom 2         R12           R10\n");
      fprintf(fp,"------------------------------------------\n");
      for(i=0; i<gNumAtomParams; i++)
      {
         for(j=0; j<gNumAtomParams; j++)
         {
            if(gHBondParams[i][j].AtomD[0])
            {
               fprintf(fp,"%4s      %4s   %11.3f   %11.3f\n",
                       gHBondParams[i][j].AtomD,
                       gHBondParams[i][j].AtomA,
                       gHBondParams[i][j].r12,
                       gHBondParams[i][j].r10);
            }
         }
      }
      fprintf(fp,"\n");
   }
}





