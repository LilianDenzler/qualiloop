/***********************************************************************
 *      Name: PimmsTypes.h                                             *
 *  Function: Declares some types                                      *
 * Copyright: (C) Oxford Molecular Limited                             *
 *---------------------------------------------------------------------*
 *    Author: J.P.Holland, Oxford Molecular Ltd.                       *
 *      Date: 22-Jul-1991                                              *
 *---------------------------------------------------------------------*
 * Modification Record                                                 *
 * Date     Inits   Comments                                           *
 ***********************************************************************/

#ifndef __PIMMSTYPES__
#define __PIMMSTYPES__

/* Define some sizes for the atom coordinates, charges etc */
#ifndef AtomCoordType
 typedef float AtomCoordType ;
#endif

#ifndef AtomChargeType
 typedef float AtomChargeType ;
#endif

#ifndef AtomPtrType
 typedef short AtomPtrType ;
#endif

#endif
