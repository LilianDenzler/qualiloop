/****************************************************************************** 
 *	    Name: ActiveMols.h                                                *
 *  Function: Molecules active ( for display, energy, writing etc ).          *
 *            These definitions were previously made in several headers.      *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: J.P.Holland, Oxford Molecular Ltd.                              *
 *      Date: 23-Aug-1990                                                     *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 ******************************************************************************/

#ifndef __ACTIVEMOLS__
#define __ACTIVEMOLS__

#define NO_MOLECULES			0
#define ALL_ACTIVE			1
#define ALL_VISIBLE			2
#define ALL_MOLECULES			3
#define ONE_MOLECULE 			4
#define CLICKED_SELECTION               5

/* -- Definitions for molecule/residue selection by clicking on an atom -- */

#define MAX_SELECT_MOLS  10
#define MAX_SELECT_RES   10
#define INITIALISE_VALUE -1

/* inactive application identifier */
#define NOTACTIVE			-1
#define NO_OPTION			-1
#define NO_ATOM				-1

#endif
