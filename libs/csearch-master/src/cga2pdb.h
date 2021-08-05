#ifndef __CGA2PDB_H__
#define __CGA2PDB_H__

/*****************************************************************************
 *      Name: cga2pdb.h                                                      *
 *  Function: Header declarations for the CGA to PDB converter               *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1992                                 *
 *---------------------------------------------------------------------------*
 *    Author: Dr. D. Walker, Oxford Molecular Ltd                            *
 *      Date: 31/03/92                                                       *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 * 09/09/92   DW         Add TRUE/FALSE and clgAtomList structure            *
 *****************************************************************************/

/* TRUE and FALSE definitions */
#ifndef FALSE
#define FALSE 0
#define TRUE  !FALSE
#endif

/* Exit codes */

#define EXITOK      0
#define BADARGS     1
#define OPENCGA     2
#define READHDR     3
#define ALLOCALT    4
#define READDUM     5
#define BADBLOCK    6
#define ALLOCX      7
#define ALLOCY      8
#define ALLOCZ      9
#define SKIPBLOCK   10
#define READBLOCK   11
#define OPENPDB     12
#define OPENTEMPPDB 13
#define RENAME      14
#define ALLOCCLG    15
#define OPENCLG     16
#define OPENTEMPCLG 17

/* General buffer space size declaration */
#define BUFSIZE     512

/* Temporary rewritten PDB file name */
#define REWRITEPDB  "cga2pdb.reftmp"
#define REWRITECLG  "cga2pdb.clgtmp"

/* Linked list of altered residue numbers and atom names */

#define ATNMLEN 6

typedef struct clgAtomList clgAtomList;
struct clgAtomList
{
   clgAtomList *last;
   clgAtomList *next;

   int  listPos;
   int  residueNumber;
   char atomName[ATNMLEN];
};

#endif
