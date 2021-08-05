/****************************************************************************** 
 *      Name: AtomTypes.h                                                     *
 *  Function: Header file - contains definitions of data types used by some   *
 *            of the functions that manipulate individual atoms - mainly      *
 *            functions in the Atoms segment.                                 *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Tessella Support Services Ltd.                      *
 *      Date: 11/04/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 ******************************************************************************/

#ifndef __ATOMTYPES__
#define __ATOMTYPES__

#define MAX_HIGHLIGHT_ATOMS		   50

/* the HighlightedAtomList type is a structure which */
/* contains a simple list of all highlighted atoms  */
typedef struct {
   AtomPtrType NHighlightedAtoms ;	      /* number of highlighted atoms */
   AtomPtrType HighlightedAtom [ MAX_HIGHLIGHT_ATOMS ] ;/* highlighted atoms */
		} HighlightedAtomList ;

#endif
