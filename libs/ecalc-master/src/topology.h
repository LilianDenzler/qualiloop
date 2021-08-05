/*************************************************************************

   Program:    ECalc
   File:       topology.h
   
   Version:    V1.5.1
   Date:       07.01.21
   Function:   Define topology structure types, etc.
   
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
   V1.4   18.05.95 Modified atom structure for shake support
                   Also added Mobile flag
   V1.5   06.02.03 Skipped
   V1.5.1 07.01.21 Skipped

*************************************************************************/
#ifndef _TOPOLOGY_H
#define _TOPOLOGY_H

/************************************************************************/
/* Includes
*/
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXEXCL  6
#define MAXCHAIN 6

/************************************************************************/
/* Structure definitions
*/
typedef struct _atom
{
   struct _atom *excl[MAXEXCL];
   REAL mass,
        pol,
        NEff,
        vdwr,
        charge,
        x,    y,    z,
        Refx, Refy, Refz;
   int  key,
        nexcl;
   BOOL Skip,
        SkipPrev,
        Mobile;
   char atom[8],
        type[8];
}  ATOM;
   
typedef struct _bond
{
   struct _bond *next;
   REAL force,
        OptLen;
   ATOM *atom1,
        *atom2;
}  BOND;

typedef struct _angle
{
   struct _angle *next;
   REAL  force,
         OptAng;
   ATOM  *atom1,
         *atom2,
         *atom3;
}  ANGLE;

typedef struct _torsion
{
   struct _torsion *next;
   REAL    force,
           OptTor;
   int     period;
   ATOM    *atom1,
           *atom2,
           *atom3,
           *atom4;
}  TORSION;

typedef struct _improper
{
   struct _improper *next;
   REAL     force,
            OptTor;
   ATOM     *atom1,
            *atom2,
            *atom3,
            *atom4;
}  IMPROPER;

typedef struct _acceptor
{
   struct _acceptor *next;
   ATOM  *atom;
}  ACCEPTOR;

typedef struct _donor
{
   struct _donor *next;
   ATOM     *atom1,
            *atom2,
            *atom3,
            *atom4;
}  DONOR;

typedef struct
{
   ATOM     **atoms;
   BOND     *bonds;
   ANGLE    *angles;
   TORSION  *torsions;
   IMPROPER *impropers;
   DONOR    *donors;
   ACCEPTOR *acceptors;
   int      *ResStart;
   char     *sequence;
   int      NRes,
            NAtoms,
            NBonds,
            NAngles,
            NTorsions,
            NImpropers,
            NDonors,
            NAcceptors;
   BOOL     Disulphide;
   char     ChainName;
}  TOPOLOGY;

typedef struct
{
   TOPOLOGY *topol[MAXCHAIN];
   int      NChains,
            NAtoms,
            NDonors;
}  MOLECULE;

typedef struct
{
   ATOM *atom,
        *antecedant;
   ATOM **neighbours;
   int  NNeighb,
        NeighbMax;
}  GRID;

/************************************************************************/
#endif
