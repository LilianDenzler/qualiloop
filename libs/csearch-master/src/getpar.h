/***********************************************************************
 *      NAME: GETPAR                                                   *
 *  FUNCTION: To point to parameters for bond etc energy calculations  *
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
 
   Variable name   Array bounds             Description
   -------------   ------------             -----------
 
       IBNDP        MAXBND      Location of bond parameters 
       IANGP        MAXANG      Location of angle parameters
       ITORP        MAXTOR      Location of proper torsion parameters
       IIMPP        MAXIMP      Location of improper torsion parameters
       IHBP         MAXHB       Location of hydrogen-bond parameters
*/
#ifndef __GETPAR_H__
#define __GETPAR_H__

extern struct {
       short ibndp[maxbnd], iangp[maxang], itorp[maxtor], 
             iimpp[maximp], ihbp[maxhb];
       } getpar;

#endif
