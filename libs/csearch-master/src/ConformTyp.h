/****************************************************************************** 
 *      Name: ConformTyp                                                      *
 *  Function: Conformations types header file - contains definition of the    *
 *            data types used by the Conformations functions.                 *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Tessella Support Services plc                       *
 *      Date: 21/05/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 19/10/90 drm     Modified for PIMMS                                        *
 ******************************************************************************/
#ifndef __CONFORMTYP__
#define __CONFORMTYP__
 
#ifndef __CONFORM__
#include "Conform.h"
#endif

/* The Conformations Bond record contains the layout of the controls and */
/* the current state of any selections made during the Bond setup dialog */
typedef struct 
	{
	short	nselectedTorsions;	/* number of atoms selected */
	short	selectedTorsions[MAX_CONF_BONDS][TORSION];	/* selected torsions */
	float	startAngle[MAX_CONF_BONDS];	/* torsion start angles */
	short	molecule;			/* molecule being scanned */
	short	ncurrent;		/* number of currently selected atoms */
	short	currentTorsion[TORSION];	/* currently selected atom pair */
	short	BSoption;		/* Bond setup option selected */
	} ConfBondRecord;

/* The Conformations Parameters record contains the    */
/* current value of the conformational search settings */
typedef struct 
	{
	float	angleIncrement;		/* torsion angle increment */
	Logical	drawEachConformer;	/* draw each conformer flag */
	Logical	energyScan;		/* energy scan flag */
	Logical	geometryScan;		/* geometry scan flag */
	short	totalconformers;	/* number of conformers in view */
	} ConfParameters;

/* The Conformations Scan record contains the layout of the controls */
/* and the current state of any selections made during the dialog    */
typedef struct 
	{
	short	nselected;		/* number of atoms so far selected */
	short	selectedAtoms[TORSION]; /* list of selected atoms */
	short	ndistances;		/* number of distances defined */
	short	distance[MAX_CONF_SCAN_VARIABLES][BOND];	/* list of atoms defining monitor distances */
	float	mindistance[MAX_CONF_SCAN_VARIABLES];		/* list of distance minima */
	float	maxdistance[MAX_CONF_SCAN_VARIABLES];		/* list of distance maxima */
	short	nangles;		/* number of angles defined */
	short	angle[MAX_CONF_SCAN_VARIABLES][ANGLE];		/* list of atoms defining monitor angles */
	float	minangle[MAX_CONF_SCAN_VARIABLES];			/* list of angle minima */
	float	maxangle[MAX_CONF_SCAN_VARIABLES];			/* list of angle maxima */
	short	ntorsions;		/* number of torsions defined */
	short	torsion[MAX_CONF_SCAN_VARIABLES][TORSION];	/* list of atoms defining monitor torsion */
	float	mintorsion[MAX_CONF_SCAN_VARIABLES];		/* list of torsion minima */
	float	maxtorsion[MAX_CONF_SCAN_VARIABLES];		/* list of torsion maxima */
	} ConfScanRecord;

/* The Conformations View record contains the layout of the controls and */
/* the current state of any selections made during the View Scan dialog */
typedef struct 
	{
	short	VSoption;	/* View scan option selected */
	short	nconformers;	/* number of conformers in view */
	short	selectedconformer;	/* selected conformer */
	char	conformers[MAX_CONFORMATIONS][5];	/* conformers list */
	} ConfViewRecord;

/* The Conformations Plot record contains the lists   */
/* of energies and angles to generate the Energy plot */
typedef struct 
	{
	Logical active;			/* is plot active */
	short	nconformers;		/* number of conformers */
	short	selectedconformer;	/* the one selected */
	float	energies[MAX_ENERGY_POINTS];	/* energy list */
	float	angles[MAX_ENERGY_POINTS];		/* angle list */
	float	maxenergy;			/* energy maximum */
	float	minenergy;			/* energy minimum */
	} ConfPlotRecord;

/* The Conformations Contour record contains the lists of  */
/* energies and angles to generate the Energy Contour plot */
typedef struct 
	{
	Logical	active;			/* active flag */
	short	nconformers[2];		/* number of conformers */
	short	selectedconformer;	/* selected conformer */
	short	selectedpoint[2];	/* selected conformer point */
	float	maxenergy;		/* energy maximum */
	float	minenergy;		/* energy minimum */
	float	levels[MAX_CONTOUR_LEVELS];	/* energy contour levels */
	short	colours[MAX_CONTOUR_LEVELS];	/* contour colour indices */
	} ConfContourRecord;

/* the ConfContourCommon type is a structure which */
/* matches the structure of the common block /CNFCTR/ */
typedef struct 
	{ float step;				/* grid step size */
	  float x[MAX_CONTOUR_POINTS + 1];	/* list of grid x-coordinates */
	  float y[MAX_CONTOUR_POINTS + 1];	/* list of grid y-coordinates */
	  float surface[MAX_CONTOUR_POINTS+1][MAX_CONTOUR_POINTS+1];
	} ConfContourCommon, *ConfContourCommonPtr;

#endif
