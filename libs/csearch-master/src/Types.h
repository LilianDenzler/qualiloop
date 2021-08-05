/***********************************************************************
 *      Name: Types.h                                                  *
 *  Function: Declares some types                                      *
 * Copyright: (C) Oxford Molecular Limited                             *
 *---------------------------------------------------------------------*
 *    Author: G C Calvert, Tessella Support Services plc               *
 *      Date: 16/07/90                                                 *
 *---------------------------------------------------------------------*
 * Modification Record                                                 *
 * Date     Inits   Comments                                           *
 * 20/08/90 JPH     Add some atom data size definitions.               *
 *                  Defining these sizes here globally will make       *
 *                  recompilation with a bigger size much easier.      *
 ***********************************************************************/
#ifndef __TYPES__
#define __TYPES__

/* define the type logical */
#ifndef Logical
 typedef unsigned char Logical ;
#endif

/* Define FORTRAN types */
typedef int    FINT4 ;
typedef float  FREAL4 ;
typedef double FREAL8 ;

#define CHAR_PER_FINT4   4
#define CHAR_PER_FREAL8  8
#define FINT4_PER_FREAL8 2

#define NULLCHAR  '\0'
#define NULLFUNC  ((void *) 0)

#endif
