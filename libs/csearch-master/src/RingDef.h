/******************************************************************************
 *      Name: RingDef.h                                                       *
 *  Function: Declarations of ring definitions and structures.                *
 * Copyright: (C) Oxford Molecular Limited                                    *
 *----------------------------------------------------------------------------*
 *    Author: J.P.Holland, Oxford Molecular Ltd.                              *
 *      Date:  4-Feb-1991                                                     *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 ******************************************************************************/

#ifndef __RINGDEF__
#define __RINGDEF__

typedef struct ringdat * RingDataPtr ;
typedef struct ringdat {	/* Structure for linked list of ring data */
	AtomPtrType   Molecule ;
	AtomPtrType   RingSize ;
	AtomPtrType * RingAtoms ;
	RingDataPtr   NextRing ;
	RingDataPtr   PrevRing ;
	RingDataPtr   LastRing ;
		} RingData ;

#endif
