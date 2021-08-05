/***********************************************************************
 *      NAME: HBONDS                                                   *
 *  FUNCTION: To declare hydrogen bonding parameters                   *
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
   Variable   Array Bounds              Description
   --------   ------------              -----------
   
   HBCUT          -          The hydrogen bond distance limit
   HBACUT         -          The hydrogen bond angle limit
   HBDON        MAXHB        The heavy atom donor for each hydrogen bond
   HBACC        MAXHB        The heavy atom acceptor for each hydrogen 
*/

#ifndef __HBONDS_H__
#define __HBONDS_H__

extern struct {
       short hbdon[maxhb], hbacc[maxhb];

       float hbcut, hbacut;

       } hbonds;

#endif
