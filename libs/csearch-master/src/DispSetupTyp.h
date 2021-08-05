/****************************************************************************** 
 *      Name: DisplaySetupTypes                                               *
 *  Function: Display Setup types header file - contains definition of the    *
 *             data types used by the Nemesis display setup functions.        *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Tessella Support Services Ltd.                      *
 *      Date: 26/03/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 * 10/10/91 KYC     Added Patom to  BondDraw structure                        *
 * 05/05/93 KJW     Make sure that DispSetup.h is always included (CR 1792)   *
 ******************************************************************************/

#ifndef __DISPSETUPTYP__
#define __DISPSETUPTYP__

#if !defined(__DISPLAY__)
#include "DispSetup.h"
#endif

/* The Display Settings record contains the current state of */
/* any settings set bythe user from the display setup menu  */
typedef struct 
	{
	Logical		Hydrogens;	/* flags if hydrogens are to be displayed */
	short		HydrogenClass;	/* hydrogen display target class */
	Logical		Dummies;	/* flags if dummies are to be displayed */
	Logical		DummyLabels;	/* flags if dummy labels are on/off */
	Logical		Labels;		/* flags if hydrogens are to be displayed */
	short		LabelClass;	/* label display target class */
	short		LabelType;	/* label type */
	short		MColourClass;	/* molecule colour display target class */
	short		MColourType;	/* molecule colour type */
	RGBColour	MRGBColour;	/* molecule colour */
	Logical		UserOrigin;	/* flags if user has defined rotation origin */
	Point3D		Origin;		/* the user origin */
	} DisplaySettings;

/* The Atom Colour record contains the layout of the controls and the */
/* current state of any selections made during the Atom Colour dialog */
typedef struct 	{
	short		nselectedAtoms;		/* number of atoms selected */
	short		selectedAtoms[MAX_AC_ATOMS];	/* selected atom to be coloured */
	RGBColour	AtomRGBColour;		/* selected atom colour */
	short		ACoption;	/* Atom Colour option selected */
	} AtomColour;

/* The Atom Type Colours record contains the current colour settings */
/* for each atom type set by the user from the display setup menu    */
typedef struct 
	{
	RGBColour	Hydrogen;		/* colour of hydrogen atoms */
	RGBColour	Carbon;			/* colour of carbon atoms */
	RGBColour	Nitrogen;		/* colour of nitrogen atoms */
	RGBColour	Oxygen;			/* colourof oxygen atoms */
	RGBColour	Phosphorus;		/* colour of phosphorus atoms */
	RGBColour	Sulphur;		/* colour of sulphur atoms */
	RGBColour	Fluorine;		/* colour of fluorine atoms */
	RGBColour	Chlorine;		/* colour of chlorine atoms */
	RGBColour	Bromine;		/* colour of bromine atoms */
	RGBColour	Iodine;			/* colour of iodine atoms */
	RGBColour	Other;			/* colour of other atoms */
	char		OtherName[3];		/* type of other atom */
	} AtomTypeColours;

/* The Origin record contains the layout of the controls and the     */
/* current state of any selections made during the Set Origin dialog */
typedef struct
        {
	short	selectedAtom;		/* selected atom - to use as origin */
	short	SOoption;		/* Set Origin option selected */
	short	SOtype;			/* Set Origin option selected */
	} OriginRecord;

/* The View Bond record contains the layout of the controls and the     */
/* current state of any selections made during the View Bond dialog */
typedef struct 
	{
	short	nAtoms;			/* number of selected atoms */
	short	newAtom;		/* first atom in bond */
	short	BondAtom[MAX_DS_VB_ATOMS]; /* atoms defining the bond */
	short	VBoption;		/* Set Origin option selected */
	} DSViewBond;

/* Bond Draw is used in the Setbdr routine */
typedef struct
	{
        Logical HighAtom ;
        Logical Padding[3];
        int RibAtom      ;
	int Molecule     ; 	
        int NumBonds     ;
        short Patom      ;  /* parent atom of bonds */
        short Padding2   ;
	AtomPtrType Bonds[MAXCON] ;
	} BondDraw ;
#endif
