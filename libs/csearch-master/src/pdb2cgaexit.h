#ifndef _PDB2CGA_EXITS_H_
#define _PDB2CGA_EXITS_H_

/*****************************************************************************
 *      Name: pdb2cgaexit.h                                                  *
 *  Function: Define exit codes for the PDB to CGA converter                 *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1992                                 *
 *---------------------------------------------------------------------------*
 *    Author: Dr D Walker, Oxford Molecular Ltd                              *
 *      Date: 17/03/92                                                       *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/

#define EXITOK      0
#define BADUSAGE    1
#define BADLOOP     2
#define BADNUMPDB   3
#define OPENCGAOUT  4
#define OPENCGAIN   5
#define EOFONCGAIN  6
#define READCGAIN   6
#define ALLOCALTER  7
#define OPENPDB     8
#define EOFONPDB    9

#define BUFSIZE     512

#ifndef FALSE
#define FALSE       0
#define TRUE        !FALSE
#endif

#define NUMPERLINE  13

#define ENERGYCOMMENT "REMARK Energy="

#endif
