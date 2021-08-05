#ifndef __OMLIB_PRO__
#define __OMLIB_PRO__

/*****************************************************************************
 *      Name: OMLIB.pro                                                      *
 *  Function: Prototype definitions for OMLIB "C" functions                  *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author: Dr. D. Walker, Oxford Molecular Ltd                            *
 *      Date: 18/01/93                                                       *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 * 01/11/93   RKW        Added prototypes for IsInteger + IsRealNumber       *
 *****************************************************************************/

/* CM */
int CIsInt(char *ctest);
void Cmljust(char *text);

Boolean IsInteger( char *text );
Boolean IsRealNumber( char *text );

#endif
