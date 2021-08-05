/****************************************************************************** 
 *      Name: Conform.h                                                       *
 *  Function: Header file - contains definitions of constants for use in the  *
 *            Conformations windows, dialogs and functions.                   *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Tessella Support Services plc                       *
 *      Date: 06/06/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 19/10/90 drm     Modified for PIMMS                                        *
 * 27/11/90 jph     Removed definitions of ATOM,BOND etc - already defined    *
 * 26/03/91 GCC     Swapped conformation save/cancel codes.                   *
 ******************************************************************************/
#ifndef __CONFORM__
#define __CONFORM__
	
/* Conformations item identifiers */
#define CONF_BOND_SETUP			0
#define CONF_GEOMETRY_SCAN		1
#define CONF_SEARCH			2
#define CONF_VIEW_SCAN			3
#define CONF_FINAL_SCAN			4
#define CONF_OPEN_SCAN			5
#define CONF_SAVE_SCAN			6
#define CONF_FINISH			7
#define CONF_MAX_OPTIONS	CONF_FINISH+1

/* Bond Setup window identifiers */
#define MAX_CONF_BONDS		  5	/* maximum number of selected bonds */
#define MAX_CONF_BOND_BUTTONS	  4	/* number of buttons on Conformations Bond setup dialog */
#define CONF_BOND_ADD		  0	/* add bond */
#define CONF_BOND_LIST		  1	/* list currently defined bonds */
#define CONF_BOND_CLEAR		  2	/* clear currently defined bonds */
#define CONF_BOND_FINISH	  3	/* end Bond Setup dialog */

/* Search setup dialog item identifiers */
#define MAX_CONF_INCREMENT 	   9	/* maximum number angle increments */
#define CONF_SEARCH_ENTER	1
#define CONF_SEARCH_CANCEL	0

/* Scan Variables window identifiers */
#define MAX_CONFORMATIONS	   1296	/* maximum number of conformations */
#define MAX_ENERGY_POINTS	    360	/* maximum number of energy plot points */
#define MAX_CONTOUR_POINTS	     36	/* maximum number of energy contour points */
#define MAX_CONTOUR_LEVELS	      5	/* maximum number of energy contour levels */
#define MAX_CONF_SCAN_VARIABLES	  5	/* maximum number of scan variables */
#define MAX_CONF_SCAN_BUTTONS	  7	/* number of buttons on Scan Variables dialog */

#define CONF_SCAN_ADD_DISTANCE		  0	/* add scan distance */
#define CONF_SCAN_ADD_ANGLE		  1	/* add scan angle */
#define CONF_SCAN_ADD_TORSION		  2	/* add scan torsion */
#define CONF_SCAN_CLEAR_DISTANCE	  3	/* clear scan distance */
#define CONF_SCAN_CLEAR_ANGLE		  4	/* clear scan angle */
#define CONF_SCAN_CLEAR_TORSION		  5	/* clear scan torsion */
#define CONF_SCAN_FINISH	          6	/* end Scan Variable Setup dialog */

#define CONF_GET_DIST_ENTER		1
#define CONF_GET_DIST_CANCEL		0

#define CONF_GET_ANGL_ENTER		1
#define CONF_GET_ANGL_CANCEL		0

/* View Scan window identifiers */
#define MAX_VIEW_CONFORMERS	  9	/* maximum number of conformers in-view */
#define MAX_CONF_VIEW_BUTTONS	  9	/* number of buttons on Scan Variables dialog */
#define CONF_VIEW_CONFORMER	  0	/* add scan distance */
#define CONF_VIEW_SATISFIED	  1	/* clear scan angle */
#define CONF_VIEW_ALL_SATISFIED	  2	/* clear scan angle */
#define CONF_VIEW_MINIMUM	  3	/* addscan angle */
#define CONF_VIEW_MAXIMUM	  4	/* add scan torsion */
#define CONF_VIEW_RANGE		  5	/* clear scan distance */
#define CONF_VIEW_PLOT		  6	/* clear scan torsion */
#define CONF_VIEW_CONTOUR_PLOT	  7	/* clear scan torsion */
#define CONF_VIEW_FINISH	  8	/* end View Scan Setup dialog */

#define CONF_VIEW_CONF_ENTER	1
#define CONF_VIEW_CONF_CANCEL   0

#define CONF_GET_ENERGY_ENTER	1
#define CONF_GET_ENERGY_CANCEL  0

#define MAX_CONF_NAME_LENGTH	5

#define CONF_CONFFILE_SAVE	1
#define CONF_CONFFILE_CANCEL	0

/* OMUPD KYC 10/12/91 add constant for Scale to zero */
/* OMUPD KYC 11/12/91 add constant for reset */
#define CONF_PLOT_RESET         3  
#define CONF_PLOT_SCALE         2
#define CONF_PLOT_ENTER		1
#define CONF_PLOT_CANCEL	0

/* Contour level dialog constants */
#define CONF_CONTOUR_LEVEL_ENTER	1
#define CONF_CONTOUR_LEVEL_CANCEL	0

#define CL_BUTTON_LEFT		100
#define CL_BUTTON_TOP		400
#define CL_BUTTON_RIGHT		130
#define CL_BUTTON_BOTTOM	350
#define CL_SHIFT		 50

#endif
