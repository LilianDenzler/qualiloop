/*************************************************************************

   Program:    ECalc
   File:       rtop.h
   
   Version:    V1.5.1
   Date:       07.01.21
   Function:   Define types for residue topology and global storage
   
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
   V1.0   30.09.94 First release version
   V1.1   12.10.94 Changes to energy.c
   V1.2   10.11.94 Bug fixes in energy.c
   V1.3   28.11.94 Bug fixes in energy.c
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Skipped
   V1.5.1 07.01.21 Skipped

*************************************************************************/
#ifndef _RTOP_H
#define _RTOP_H

/************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXRESTYPE 30
#define MAXEXCL     6


/************************************************************************/
/* Structure definitions
*/
typedef struct _atomtop
{
   struct _atomtop *next;
   REAL            charge;
   int             nexcl;
   char            atom[8],
                   type[8],
                   excl[MAXEXCL][8];
}  ATOMTOP;

typedef struct _bondtop
{
   struct _bondtop *next;
   char            atom1[8],
                   atom2[8];
}  BONDTOP;

typedef struct _angletop
{
   struct _angletop *next;
   char             atom1[8],
                    atom2[8],
                    atom3[8];
}  ANGLETOP;

typedef struct _torsiontop
{
   struct _torsiontop *next;
   char               atom1[8],
                      atom2[8],
                      atom3[8],
                      atom4[8];
}  TORSIONTOP;

typedef struct _impropertop
{
   struct _impropertop *next;
   char                atom1[8],
                       atom2[8],
                       atom3[8],
                       atom4[8];
}  IMPROPERTOP;

typedef struct _donortop
{
   struct _donortop *next;
   char             atom1[8],
                    atom2[8],
                    atom3[8],
                    atom4[8];
}  DONORTOP;

typedef struct _acceptortop
{
   struct _acceptortop *next;
   char                atom[8];
}  ACCEPTORTOP;

typedef struct
{
   char        name[8];
   REAL        charge;
   ATOMTOP     *AtomTop;
   BONDTOP     *BondTop;
   ANGLETOP    *AngleTop;
   TORSIONTOP  *TorsionTop;
   IMPROPERTOP *ImproperTop;
   DONORTOP    *DonorTop;
   ACCEPTORTOP *AcceptorTop;
}  RESIDUETOP;


/************************************************************************/
/* Globals
*/
#ifdef READRTOP_MAIN    /*--------------- Define globals ---------------*/
RESIDUETOP              gRTop[MAXRESTYPE];
int                     gNumResRTop = (-1);

#else                   /*------------- Reference globals --------------*/
extern RESIDUETOP       gRTop[MAXRESTYPE];
extern int              gNumResRTop;

#endif                  /*----------------- End globals ----------------*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
#endif











