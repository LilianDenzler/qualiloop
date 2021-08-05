#ifndef __OMLIB_H__
#define __OMLIB_H__

/*****************************************************************************
 *      Name: OMLIB.h                                                        *
 *  Function: General header definitions for OMLIB                           *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author: Dr. D. Walker, Oxford Molecular Ltd                            *
 *      Date: 18/01/93                                                       *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 * 01/11/93   RKW        Added NemTypes.h for new OMLIB.pro contents         *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "NemTypes.h"

/* Make sure TRUE and FALSE are defined */
#define FALSE 0
#define TRUE  1

/* Make sure the function prototypes are defined */
#include "OMLIB.pro"

#endif
