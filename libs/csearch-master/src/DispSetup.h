/****************************************************************************** 
 *	Name: DisplaySetup                                                    *
 *  Function: Display Setup header file - contains definition of the data     *
 *            items used by PIMMS in the Display Setup Application.           *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: D.R. Marsh, Tessella Support Services plc                       *
 *      Date: 08/08/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 * 15/01/91 GCC     Removed stick and ball/stick options                      *
 * 08/04/91 GCC     Added center molecule option.                             *
 * 15/07/92 CAB     Added DT_BALLSTICK for ball & stick option                *
 ******************************************************************************/

#ifndef __DISPLAY__
#define __DISPLAY__
#endif

#define HYDROGENS_OFF	0
#define HYDROGENS_ON 	1

/* Set Origin window identifiers and constants */

#define SO_USER_DEFINED		  1	/* user defined origin flag */
#define SO_DEFAULT		  0	/* default origin flag */

#define SO_SET_ORIGIN		  0	/* set origin flag */
#define SO_RESET_ORIGIN		  1	/* reset to default origin flag */
#define SO_CENTRE      		  2	/* centre molecules */
#define SO_FINISH		  3	/* end Set Origin dialog */

/* View Bond identifiers and constants */

#define VB_VIEW_BOND 		  0	/* View Bond action */
#define VB_CLEAR       		  1	/* clear highlighted atoms */
#define VB_FINISH		  2	/* end View Bond dialog */

#define MAX_DS_VB_ATOMS		2	/* max atoms that can be selected */

/* Atom Colour identifiers and constants */

#define AC_SET_SELECTION	0	/* set selected atom colour flags */
#define AC_SET_TYPES		1	/* set colours by atom type flags */
#define AC_FINISH		2	/* end Atom Colour dialog */

#define MAX_AC_ATOMS		20
#define MAX_AC_BUTTONS		  3	/* number of buttons on atom colour dialog */

/* Atom colour dialog item identifiers */
#define iDACHydrogen		  10	/* hydrogen colour user item */
#define iDACCarbon		   9	/* carbon colour user item */
#define iDACNitrogen		   8	/* nitrogen colour user item */
#define iDACOxygen		   7	/* oxygen colour user item */
#define iDACPhosphorus		   6	/* phosphorus colour user item */
#define iDACSulphur		   5	/* sulphur colour user item */
#define iDACFluorine		   4    /* fluorine colour user item */
#define iDACChlorine		   3	/* chlorine colour user item */
#define iDACBromine		   2	/* bromine colour user item */
#define iDACIodine		   1	/* iodine colour user item */
#define iDACOther		   0	/* other element colour user item */
#define iDACOtherName		  25	/* other element name item */

#define AC_BUTTON_LEFT		180	/* left corner of first button */
#define AC_BUTTON_TOP		105	/* top corner of first button */
#define AC_BUTTON_RIGHT		240	/* right corner of first button */
#define AC_BUTTON_BOTTOM	 85	/* bottom corner of first button */
#define AC_SHIFT		 40	/* inter-button distance */
 
/* Display type constants */

#define	DT_SPACEFILL		0	/* Spacefill display */
#define	DT_BALLSTICK		1	/* Ball & Stick display */
#define	DT_POTENTIAL		2	/* Potential display */
#define DT_FASTFILL             3       /* Fast space "fill" display */
#define DT_STICK                4       /* Stick Model display */
#define DT_FINISH		5	/* end display type dialog flag */

