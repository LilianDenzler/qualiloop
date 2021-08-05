/****************************************************************************** 
 *      Name: ResidueTypes.h                                                  *
 *                                                                            *
 *  Function: Header file - contains definitions of data types used by some   *
 *            of the functions that manipulate residues.                      *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: D.R. Marsh, Tessella Support Services plc                       *
 *      Date: 31/12/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 ******************************************************************************/

#ifndef __RESIDUETYPES__
#define __RESIDUETYPES__

#define MAX_HIGHLIGHT_RESIDUES  2

/* the HighlightedResidueList type is a structure which */
/* contains a simple list of all highlighted Residues  */
typedef struct 
	{ short NHighlightedResidues;	/* number of highlighted residues */
	  short HighlightedResidue[MAX_HIGHLIGHT_RESIDUES];	/* list of highlighted residues */
	} HighlightedResidueList;

#endif
