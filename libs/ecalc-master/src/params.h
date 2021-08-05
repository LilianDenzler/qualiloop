/*************************************************************************

   Program:    ECalc
   File:       params.h
   
   Version:    V1.5.1
   Date:       07.01.21
   Function:   Variable type definitions and global storage for parameter
               data
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
   V0.1   25.08.94 Original
   V1.0   30.09.94 First release version
   V1.1   11.10.94 Added statistical residue pseudo-energy
   V1.2   10.11.94 Bug fixes in energy.c
   V1.3   28.11.94 Bug fixes in energy.c
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Skipped
   V1.5.1 07.01.21 Skipped

*************************************************************************/
#ifndef _PARAMS_H
#define _PARAMS_H

/************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBONDP        100
#define MAXANGLEP       300
#define MAXTORSIONP      75
#define MAXIMPROPERP     75
#define MAXHBONDP       200
#define MAXATOMP         75
#define MAXNONBONDP     (MAXATOMP * MAXATOMP)

/************************************************************************/
/* Structure definitions
*/
typedef struct 
{
   REAL force,
        OptLen;
   char atom1[8],
        atom2[8];
}  BONDPARAM;

typedef struct 
{
   REAL force,
        OptAng;
   char atom1[8],
        atom2[8],
        atom3[8];
}  ANGLEPARAM;

typedef struct 
{
   REAL force,
        period,
        OptTor;
   char atom1[8],
        atom2[8],
        atom3[8],
        atom4[8];
}  TORSIONPARAM;

typedef struct 
{
   REAL force,
        OptTor;
   char atom1[8],
        atom2[8],
        atom3[8],
        atom4[8];
}  IMPROPERPARAM;

typedef struct 
{
   REAL r12,
        r10;
   char AtomD[8],
        AtomA[8];
}  HBONDPARAM;

typedef struct
{
   REAL r6,
        r12;
   char atom1[8],
        atom2[8];
}  NONBONDPARAM;

typedef struct
{
   REAL mass,
        pol,
        NEff,
        vdwr;
   int  key;
   char atom[8];
}  ATOMPARAM;

typedef struct _resparam
{
   struct _resparam *next;
   REAL             *values;
   REAL             maxval,
                    scale;
   int              nbox;
   char             distal[8],
                    restype;
}  RESPARAM;

/************************************************************************/
/* Globals
*/
#ifdef READPARAMS_MAIN  /*--------------- Define globals ---------------*/
BONDPARAM      gBondParams[MAXBONDP];
ANGLEPARAM     gAngleParams[MAXANGLEP];
TORSIONPARAM   gTorsionParams[MAXTORSIONP];
IMPROPERPARAM  gImproperParams[MAXIMPROPERP];
NONBONDPARAM   **gNonBondParams = NULL;
HBONDPARAM     **gHBondParams   = NULL;
ATOMPARAM      gAtomParams[MAXATOMP];
RESPARAM       *gResParam = NULL;

int            gNumBondParams     = 0,
               gNumAngleParams    = 0,
               gNumTorsionParams  = 0,
               gNumImproperParams = 0,
               gNumAtomParams     = 0,
               gNumNonBondParams  = 0,
               gNumHBondParams    = 0;

#else                   /*------------- Reference globals --------------*/
extern BONDPARAM      gBondParams[MAXBONDP];
extern ANGLEPARAM     gAngleParams[MAXANGLEP];
extern TORSIONPARAM   gTorsionParams[MAXTORSIONP];
extern IMPROPERPARAM  gImproperParams[MAXIMPROPERP];
extern NONBONDPARAM   **gNonBondParams;
extern HBONDPARAM     **gHBondParams;
extern ATOMPARAM      gAtomParams[MAXATOMP];
extern RESPARAM       *gResParam;

extern int            gNumBondParams,
                      gNumAngleParams,
                      gNumTorsionParams,
                      gNumImproperParams,
                      gNumAtomParams,
                      gNumNonBondParams,
                      gNumHBondParams;

#endif                  /*----------------- End globals ----------------*/

/************************************************************************/
/* Prototypes
*/
#include "GetParam.p"

/************************************************************************/
#endif
