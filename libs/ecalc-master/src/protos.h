/*************************************************************************

   Program:    ECalc
   File:       protos.h
   
   Version:    V1.5.1
   Date:       07.01.21
   Function:   Include prototype files
   
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
   V1.0   30.09.94 First release version
   V1.1   12.10.94 Changes to energy.c
   V1.2   10.11.94 Bug fixes in energy.c
   V1.3   28.11.94 Bug fixes in energy.c
   V1.4   18.05.95 Added shake.p
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Skipped
   V1.5.1 07.01.21 Skipped

*************************************************************************/
/* Includes
*/
#ifndef _PROTOS_H
#define _PROTOS_H

void PrintParams(FILE *fp);
BOOL ReadParams(char *filename);
BOOL ReadRTop(char *filename);
void PrintRTop(FILE *out);
void FreeRTop(void);
int  FindRTop(char *resnam);

#include "BuildTop.p"
#include "GetParam.p"
#include "energy.p"
#include "ReadStructure.p"
#include "main.p"
#include "shake.p"

#endif
