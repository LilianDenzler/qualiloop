/*************************************************************************

   Program:    ECalc
   File:       GetParam.c
   
   Version:    V1.5.1
   Date:       07.01.21
   Function:   Get offsets into the parameter arrays.
   
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
   V0.1   30.08.94 Original
   V1.0   30.09.94 First release version
   V1.1   12.10.94 Changes to energy.c
   V1.2   10.11.94 Bug fixes in energy.c
   V1.3   28.11.94 Bug fixes in energy.c
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Skipped
   V1.5.1 07.01.21 Skipped

*************************************************************************/
/* Includes
*/
#define GETPARAM_MAIN

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/parse.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

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
/*>int FindAtomParam(char *type)
   -----------------------------
   Search the parameter array for the offset to a given atom type

   30.08.94 Original    By: ACRM
*/
int FindAtomParam(char *type)
{
   int i;

   for(i=0; i<gNumAtomParams; i++)
   {
      if(!strncmp(gAtomParams[i].atom, type, 4))
         return(i);
   }

   return(-1);
}

/************************************************************************/
/*>int FindBondParam(char *type1, char *type2)
   -------------------------------------------
   Search the parameter array for the offset to a given bond type

   30.08.94 Original    By: ACRM
*/
int FindBondParam(char *type1, char *type2)
{
   int i;
   
   for(i=0; i<gNumBondParams; i++)
   {
      if((!strncmp(type1, gBondParams[i].atom1, 4) &&
          !strncmp(type2, gBondParams[i].atom2, 4)) ||
         (!strncmp(type1, gBondParams[i].atom2, 4) &&
          !strncmp(type2, gBondParams[i].atom1, 4)))
      {
         return(i);
      }
   }
   return(-1);
}

/************************************************************************/
/*>int FindAngleParam(char *type1, char *type2, char *type3)
   ---------------------------------------------------------
   Search the parameter array for the offset to a given angle type

   30.08.94 Original    By: ACRM
*/
int FindAngleParam(char *type1, char *type2, char *type3)
{
   int i;
   
   for(i=0; i<gNumAngleParams; i++)
   {
      if((!strncmp(type1, gAngleParams[i].atom1, 4) &&
          !strncmp(type2, gAngleParams[i].atom2, 4) &&
          !strncmp(type3, gAngleParams[i].atom3, 4)) ||
         (!strncmp(type1, gAngleParams[i].atom3, 4) &&
          !strncmp(type2, gAngleParams[i].atom2, 4) &&
          !strncmp(type3, gAngleParams[i].atom1, 4)))
      {
         return(i);
      }
   }
   return(-1);
}

/************************************************************************/
/*>int FindTorsionParam(char *type1, char *type2, char *type3, 
                        char *type4)
   -----------------------------------------------------------
   Search the parameter array for the offset to a given torsion type

   31.08.94 Original    By: ACRM
*/
int FindTorsionParam(char *type1, char *type2, char *type3, char *type4)
{
   int i;
   
   for(i=0; i<gNumTorsionParams; i++)
   {
      if(!strncmp(gTorsionParams[i].atom1, "X   ", 4))
      {
         if((!strncmp(type2, gTorsionParams[i].atom2, 4) &&
             !strncmp(type3, gTorsionParams[i].atom3, 4)) ||
            (!strncmp(type2, gTorsionParams[i].atom3, 4) &&
             !strncmp(type3, gTorsionParams[i].atom2, 4)))
         {
            return(i);
         }
      }
      else
      {
         if((!strncmp(type1, gTorsionParams[i].atom1, 4) &&
             !strncmp(type2, gTorsionParams[i].atom2, 4) &&
             !strncmp(type3, gTorsionParams[i].atom3, 4) &&
             !strncmp(type4, gTorsionParams[i].atom4, 4)) ||
            (!strncmp(type1, gTorsionParams[i].atom4, 4) &&
             !strncmp(type2, gTorsionParams[i].atom3, 4) &&
             !strncmp(type3, gTorsionParams[i].atom2, 4) &&
             !strncmp(type4, gTorsionParams[i].atom1, 4)))
         {
            return(i);
         }
      }
      
   }
   return(-1);
}

/************************************************************************/
/*>int FindImproperParam(char *type1, char *type2, char *type3, 
                         char *type4)
   ------------------------------------------------------------
   Search the parameter array for the offset to a given improper type

   31.08.94 Original    By: ACRM
*/
int FindImproperParam(char *type1, char *type2, char *type3, char *type4)
{
   int i;
   
   for(i=0; i<gNumImproperParams; i++)
   {
      if(!strncmp(gImproperParams[i].atom2, "X   ", 4))
      {
         if((!strncmp(type1, gImproperParams[i].atom1, 4) &&
             !strncmp(type4, gImproperParams[i].atom4, 4)) ||
            (!strncmp(type1, gImproperParams[i].atom4, 4) &&
             !strncmp(type4, gImproperParams[i].atom1, 4)))
         {
            return(i);
         }
      }
      else
      {
         if((!strncmp(type1, gImproperParams[i].atom1, 4) &&
             !strncmp(type2, gImproperParams[i].atom2, 4) &&
             !strncmp(type3, gImproperParams[i].atom3, 4) &&
             !strncmp(type4, gImproperParams[i].atom4, 4)) ||
            (!strncmp(type1, gImproperParams[i].atom4, 4) &&
             !strncmp(type2, gImproperParams[i].atom3, 4) &&
             !strncmp(type3, gImproperParams[i].atom2, 4) &&
             !strncmp(type4, gImproperParams[i].atom1, 4)))
         {
            return(i);
         }
      }
      
   }
   return(-1);
}

