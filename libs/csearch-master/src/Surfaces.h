/****************************************************************************** 
 *      Name: Surfaces.h                                                      *
 *  Function: Header file - contains definitions of constants for use in the  *
 *            Connolly surface windows, dialogs and functions.                *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Tessella Support Services plc                       *
 *      Date: 21/05/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 03/10/90 drm     Modified for PIMMS                                        *
 ******************************************************************************/

#ifndef __SURFACES__
#define __SURFACES__
	
/* Surfaces menu button identifiers */
#define SURF_GENERATE			0
#define SURF_EXCLUDE			1
#define SURF_DISPLAY			2
#define SURF_DELETE			3
#define SURF_OPEN			4
#define SURF_SAVE			5
#define SURF_FINISH			6
#define SURF_MAX_OPTIONS SURF_FINISH+1

/* Surfaces setup identifiers */
#define SURF_SETUP_GENERATE		1
#define SURF_SETUP_CANCEL		0

/* Surfaces display identifiers */
#define SURF_DISP_ENTER			1
#define SURF_DISP_CANCEL		0

/* Surfaces potential dialog identifiers */
#define SURF_POTENTIAL_ENTER		1
#define SURF_POTENTIAL_CANCEL		0

/* Surfaces save file dialog identifiers */
#define SURF_SAVEFILE_SAVE		1
#define SURF_SAVEFILE_CANCEL		0

/* Array and selection limits */
#define MAX_SURF_NAME			31
#define MAX_SURF_FILES			 7
#define MAX_SURF      		         2	/* maximum number of surfaces */
/* drm 26/01/91 update number of dots from 2000 to 10000 */
#define MAXDOT				10000	/* maximum number of surface points */

/* Limits for parameter settings */
#define MIN_DOT_DENSITY			  1.0
#define MAX_DOT_DENSITY			100.0
#define DEFAULT_DOT_DENSITY		  5.0
#define MIN_PROBE_RADIUS		  0.0
#define MAX_PROBE_RADIUS		  3.0
#define DEFAULT_PROBE_RADIUS		  1.4
#define MIN_VDW_EXPAND			  0.5
#define MAX_VDW_EXPAND			  2.5
#define DEFAULT_VDW_EXPAND		  1.0

/* Surface COMMON block names */
#define SURFX				 "/SURFX/"
#define SURFCS				 "/SURFCS/"
#define SURFRS				 "/SURFRS/"
#define SURFES				 "/SURFES/"

/* Surface file name */
#define SURFACE_FILE		"surface.lst"

/* Surface Colour modes */
#define MOLECULE_COLOUR			 0
#define ATOM_COLOUR			 1
#define ESPOTL_COLOUR			 2

/* Surface Types */
#define CONTACT				 0
#define REENTRANT			 1

/* inactive application identifier */
#define NOTACTIVE			-1
#define NO_OPTION			-1

/* Surface Colour Constants */
#define MAX_COLOUR_LEVELS		 5

/* Exclude Atom dialog constants */
#define SURF_EA_EXCLUDE		   0	/* exclude button ID */
#define SURF_EA_CLEAR		   1	/* clear button ID */
#define SURF_EA_FINISH		   2	/* finish button ID */
#define MAX_SURF_EA_ATOMS	  20	/* maximum number of atoms that can be selected */

/* Potential colour dialog constants */
#define PC_BUTTON_LEFT		100	/* left corner of first button */
#define PC_BUTTON_TOP		350	/* top corner of first button */
#define PC_BUTTON_RIGHT		130	/* right corner of first button */
#define PC_BUTTON_BOTTOM	300	/* bottom corner of first button */
#define PC_SHIFT		 50	/* inter-button distance */

#endif
