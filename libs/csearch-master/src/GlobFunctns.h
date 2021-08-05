/****************************************************************************** 
 *      Name: GlobFunctns.h                                                   *
 *  Function: Declaration of variables for functions accessed through global  *
 *            variables, rather than directly.
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: J.P.Holland, Oxford Molecular Ltd.                              *
 *      Date: 27-Mar-1991                                                     *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 ******************************************************************************/
#ifndef __GLOBFUNCTNS__
#define __GLOBFUNCTNS__

short ( * fProcMessage ) (
short     CancelButton ,
short     Severity ,
char    * Message ) ;

#endif
