/***********************************************************************
 *      Name: TrueFalse                                                *
 *  Function: Declares some types                                      *
 * Copyright: (C) Oxford Molecular Limited                             *
 *---------------------------------------------------------------------*
 *    Author: G C Calvert, Tessella Support Services plc               *
 *      Date: 16/07/90                                                 *
 *---------------------------------------------------------------------*
 *    Inputs:                                                          *
 *   Returns:                                                          *
 *   Globals:                                                          *
 *---------------------------------------------------------------------*
 * Modification Record                                                 *
 * Date     Inits   Comments                                           *
 * dd/mm/yy                                                            *
 ***********************************************************************/

#ifndef __TRUEFALSE__
#define __TRUEFALSE__

enum { false , true } ;

#define PIMMS_FALSE 0
#define PIMMS_TRUE 1

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#endif
