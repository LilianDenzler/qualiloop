/***********************************************************************
 *      NAME: DBG                                                      *
 *  FUNCTION: To declare the DeBuG level parameters for CONGEN         *
 * COPYRIGHT: (C) OXFORD MOLECULAR LIMITED, 1992                       *
 *---------------------------------------------------------------------*
 *    AUTHOR: Robert Williams                                          *
 *      DATE: 05/10/92                                                 *
 *---------------------------------------------------------------------*
 *    INPUTS:                                                          *
 *   OUTPUTS:                                                          *
 *    LOCALS:                                                          *
 *   GLOBALS:                                                          *
 *     CALLS:                                                          *
 *---------------------------------------------------------------------*
 * MODIFICATION RECORD                                                 *
 * DD/MM/YY   INITS   COMMENTS                                         *
 ***********************************************************************

Declarations */

#ifndef __DBG_H__
#define __DBG_H__

extern struct {
       int cgen, clschn, alloc, allstk, allhp;
       } dbg;

#endif
