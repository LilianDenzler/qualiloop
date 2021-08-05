/******************************************************************************
 *      Name: Restraints.h                                                    *
 *  Function: Declarations of restraints defintiions and structures.          *
 * Copyright: (C) Oxford Molecular Limited                                    *
 *----------------------------------------------------------------------------*
 *    Author: J.P.Holland, Oxford Molecular Ltd.                              *
 *      Date:  5-Mar-1991                                                     *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 ******************************************************************************/

#ifndef __RESTRAINTS__
#define __RESTRAINTS__

/* Define the types of restraint and constraints allowed */
/* NB The only restraints thus far supported are in eureka as follows : */
/*     fixed residues, fixed molecules and restrained torsions */
enum { RESTRAIN_START = 100 ,		/* start */
       RESTRAIN_BOND ,			/* bond length restraint */
       RESTRAIN_ANGLE ,			/* bond angle restraint */
       RESTRAIN_TORSION ,		/* torsion angle restraint */
       RESTRAIN_DISTANCE ,		/* distance restraint */

       RESTRAIN_FINISH } ;
enum { CONSTRAIN_START = 200 ,		/* start */
       CONSTRAIN_BOND ,			/* bond length constraint */
       CONSTRAIN_ANGLE ,		/* bond angle constraint */
       CONSTRAIN_TORSION ,		/* torsion angle constraint */
       CONSTRAIN_DISTANCE ,		/* distance constraint */

       CONSTRAIN_FINISH } ;
enum { FIXED_START = 300 ,		/* start */
       FIXED_ATOM ,			/* fixed atom constraint */
       FIXED_GROUP ,			/* fixed group constraint */
       FIXED_RESIDUE ,			/* fixed residue constraint */
				/* Pimms does not understand residues */
				/* so this is just another fixed group */
       FIXED_MOLECULE ,			/* fixed molecule constraint */

       FIXED_FINISH } ;


typedef struct restraints {	/* Structure for a restraint record */
        int           Type ;		/* restraint type */
	AtomPtrType * Atoms ;		/* pointer to atoms involved */
	double        Value1 ;		/* value of restraint */
	double        Value2 ;		/* weighting for restraint */
					/* or delta for constraint */
		} RestList , * RestListPtr ;

#endif
