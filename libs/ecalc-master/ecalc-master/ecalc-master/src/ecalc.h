/*************************************************************************

   Program:    ECalc
   File:       ecalc.h
   
   Version:    V1.5.2
   Date:       05.03.21
   Function:   General type definitions, defines, macros and globals
   
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
   V1.1   12.10.94 Changes to energy.c
   V1.2   10.11.94 Bug fixes in energy.c
   V1.3   28.11.94 Bug fixes in energy.c. Added global gConfNum
   V1.4   18.05.95 Shake & Relax support
   V1.5   06.02.03 Skipped
   V1.5.1 07.01.21 Skipped
   V1.5.2 05.03.21 Added gAtomError

*************************************************************************/
#ifndef _ECALC_H
#define _ECALC_H

/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "bioplib/pdb.h"
#include "bioplib/ErrStack.h"

#include "topology.h"
#include "params.h"

/************************************************************************/
/* Defines and macros
*/
#if (unix || __unix__ || msdos || __msdos__)
#define DATADIR   "ECALCDATA"     /* Environment variable name          */
#else
#define DATADIR   "ECALCDATA:"
#endif

#define PARAMFILE "params.dat"
#define TOPFILE   "restop.dat"

#define MAXBUFF 160
#define ERRBUFF 240

#define MAXRELAXITER 1            /* Max iterations for RelaxStructure()*/
#define MAXVDWITER   1            /* Max iterations for VDW Relaxation  */
#define MAXSHAKEITER 1000         /* Max iterations for SHAKE           */

#define STATE_NONE      0
#define STATE_POTENTIAL 1
#define STATE_SHOW      2

#define CHECK_STATE(x) if(state == STATE_NONE)                       \
                       {  sprintf(gError,"%s must be specified after \
POTENTIAL or DISPLAY",(x));                                          \
                          StoreError("ParseControlFile()",gError);   \
                          return(FALSE);                             \
                       }

/* VALID(ATOM *p)
   --------------
   Checks that p points to a valid atom, not a dummy

   13.09.94 Original    By: ACRM
   15.09.94 Added check on p!=NULL
*/
#define VALID(p) (((p) != NULL)          &&           \
                  ((p)->x < (REAL)9999.0 ||           \
                   (p)->y < (REAL)9999.0 ||           \
                   (p)->z < (REAL)9999.0))

typedef struct
{
   REAL GridCut,
        NonBondCut,
        eta,
        CutOnHB,
        CutOffHB,
        CutOnHBAng,
        CutOffHBAng,
        VDWTol,
        ShakeTol;
   BOOL ConstDielectric,
        DoRelax;
}  EPARAMS;

typedef struct
{
   REAL BondScale,
        AngleScale,
        TorsionScale,
        ImproperScale,
        VDWAScale,
        VDWRScale,
        HBondScale,
        ResidueScale,
        ElectScale;
   int  NCache,
        Regrid;
   BOOL PrintParams,
        PrintRTop,
        ShowBonds,
        ShowAngles,
        ShowTorsions,
        ShowImpropers,
        ShowVDWA,
        ShowVDWR,
        ShowElect,
        ShowHBonds,
        ShowResidue,
        bonds,
        angles,
        torsions,
        impropers,
        vdwa,
        vdwr,
        hbonds,
        elect,
        residue,
        Timings,
        Debug,
        AutoDisulphide;
}  FLAGS;

typedef struct _zone
{
   struct _zone *next;
   TOPOLOGY     *topol;
   int          atom1,
                atom2;
   BOOL         sc;
   char res1[16],
        res2[16];
}  ZONE;

typedef struct
{
   REAL energy;
   int  confnum;
}  CACHE;

/************************************************************************/
/* Globals
*/
#ifdef ECALC_MAIN       /*--------------- Define globals ---------------*/
char   gError[ERRBUFF];
GRID   *gGrid     = NULL,
       *gHBGrid   = NULL;
ZONE   *gUserSS   = NULL,
       *gZone     = NULL,
       *gIgnore   = NULL;
CACHE  *gCache    = NULL;
int    gConfNum   = 0;
BOOL   gAtomError = FALSE;

#else                   /*------------- Reference globals --------------*/
extern char   gError[ERRBUFF];
extern GRID   *gGrid,
              *gHBGrid;
extern ZONE   *gUserSS,
              *gZone,
              *gIgnore;
extern CACHE  *gCache;
extern int    gConfNum;
extern BOOL   gAtomError;

#endif                  /*----------------- End globals ----------------*/

/************************************************************************/
/* Prototypes
*/
#include "protos.h"

/************************************************************************/
#endif
