/****************************************************************************** 
 *      Name: Globals.h                                                       *
 *  Function: Header file to declare global variables. The "g" prefix is used *
 *            to emphasize the global nature of a variable.                   *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Tessella Support Services Ltd.                      *
 *      Date: 04/12/89                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 20/08/90 drm     Added to for PIMMS                                        *
 * 16/11/90 GCC     Added global double redraw count                          *
 * 02/12/90 drm     Added global protein settings                             *
 * 24/12/90 drm     Added new COMMON for internal residue numbers             *
 * 25/12/91 DW      Added new COMMON for the display flags                    *
 * 18/02/91 DW      Added new COMMON for the CGF key display                  *
 * 03/08/91 JPH     Remove PGL global variables to PGLGlobal.h                *
 * 18/10/91 KYC     Added global for count of visible atoms                   *
 * 27/05/92 KYC     Added global for edit labels                              *
 * 16/07/92 KYC     Added Cosmic COMMON structures                            *
 * 17/05/92 KJW     Moved COSMIC globals to cosglbls.c to make them an        *
 *                  integral part of the Cosmic library (CR 1798)             *
 ******************************************************************************/

#ifndef __GLOBALS__
#define __GLOBALS__

#ifndef __FORTRANTYPES__
#include "FORTRANTypes.h"
#endif

#ifndef __ATOMTYPES__
#include "AtomTypes.h"
#endif

#ifndef __RESIDUETYPES__
#include "ResidueTypes.h"
#endif

#ifndef __STREDITTYPES__
#include "StrEditTypes.h"
#endif

#ifndef __SURFACESTYP__
#include "SurfacesTyp.h"
#endif

#ifndef __CONFORMTYP__
#include "ConformTyp.h"
#endif

#ifndef __PROTEINTYP__
#include "ProteinTyp.h"
#endif

#ifndef __COSMICTYPES__
/* omupd kyc 16/07/92, new cosmic stuff, structures now in CosmicTYPES */
#include "COSMICTypes.h"
#endif

short   gCalcBondTarget ;         /* bond recalc.n target class */
short   gAutoAddTarget ;         /* Auto Add target class */
short   gTotalVisible;                          /* total visible atoms */
short   gSolidFill;                             /* flag for solid fill */
short   gBallAndStick;                          /* flag for Ball and stick */


/* pointers to the FORTRAN COMMON blocks */
AtomCoordCommonPtr    gAtomCoordinates;   /* pointer to Atom Coords */
AtomFlagsCommonPtr    gAtomFlags;      /* pointer to Atom Properties */
AtomPropCommonPtr     gAtomProperties;   /* pointer to Atom Properties */
AtomConnectCommonPtr  gAtomConnect;      /* pointer to Atom Connects */
AtomLabelsCommonPtr   gAtomLabels;      /* pointer to Atom Labels */
ResidueNumbersCommonPtr   gResidueNumbers;   /* pointer to internal res nos*/
BondDrawCommonPtr     gBondDraw;      /* pointer to bond draw list */
MoleculeListCommonPtr     gMoleculeList;      /* pointer to Molecule List */
MoleculeCentreCommonPtr   gMoleculeCentre;   /* pointer to Molecule Centre */
MoleculeFlagsCommonPtr    gMoleculeFlags;      /* pointer to Molecule Flags */
MoleculeTitlesCommonPtr   gMoleculeTitles;   /* pointer to Molecule Titles */
DummyCommonPtr      gDummy;         /* pointer to dummy atom defn */
PlaneCommonPtr      gPlane;         /* pointer to plane defns */
IsolatedAtomCommonPtr   gIsolatedAtoms;      /* pointer to isolated atoms */
NewFragmentCommonPtr   gNewFragment;      /* pointer to New Frag COMMON */
DummiesRecord      gDummiesStatus;      /* the dummies status record */
GeomCalcRecord      gGCStatus;      /* the geometry calcn record */
BondRotationRecord   gBRStatus;      /* bond rotation record */
FittingRecord      gFitStatus;      /* fitting record */
MonitorRecord      gMonitor;      /* monitor record */
ChargeRecord      gChargeStatus;      /* charge record */
EnergyCalcRecord   gEnergyStatus;      /* energy record */
ChiralityRecord         gChiralityStatus;   /* chirality record */
HighlightedAtomList   gHighlightedAtoms;   /* highlighted atoms record */
ConfBondRecord          gConfBondStatus;        /* Conformation Bond Record */
ConfParameters          gConfSettings;          /* Conform Search Paremeters */
ConfScanRecord          gConfScanStatus;        /* Conform Scan Variables */
ConfViewRecord          gConfViewStatus;        /* Conform Scan View record */
ConfPlotRecord          gConfPlot;              /* Conform energy plot record */
ConfContourRecord       gConfContour;           /* Conform contour plt record */
DisplaySettings      gDisplayStatus;
RGBColour           gBlack, gWhite, gRed, gGreen, gBlue, gYellow,
         gCyan, gMagenta, gGrey;
AtomColour      gAColourStatus;      /* atom colour dialog status */
AtomTypeColours      gATColours;      /* atom type colour status */
OriginRecord      gOriginStatus;      /* current origin status */
DSViewBond      gDSVBStatus;       /* view bond status */

SEResetRecord      gSEReset;
SEDeleteAtom      gSEDAStatus;
SEDeleteBond      gSEDBStatus;
SEDeleteFrag      gSEDFStatus;
SEMakeBond      gSEMBStatus;
SERetypeAtom      gSERAStatus;
SEAddHydrogens      gSEAHStatus;
SEAddFragment      gSEAFStatus;
SEChangeGroup      gSECGStatus;
SELabelAtom             gSEELStatus;            /*Label Edit status */
SETypeAtom              gSEETStatus;            /*Label Type status */

SurfaceSettings         gSurfaceStatus;
SurfExcludeAtom         gSurfEAStatus;

/* pointers to surface common blocks */
SurfacePtr              gContact;      /* contact surface common */
SurfacePtr              gReentrant;      /* reentrant surface common */
SurfaceLookupPtr        gSurfaceXref;      /* surface lookup common */
SurfacePropertyPtr      gSurfaceESPot;      /* surface property common */
SurfaceColour           gESPotColour;      /* surface property colour */

/* pointers to conformations common blocks */
ConfContourCommonPtr gConfSurfaceCtr;

/* protein structures */
SeqBuildSettings   gSeqBuildStatus;
PBAddResidue      gPBARStatus;
PBDeleteResidue      gPBDRStatus;
PBInsertResidue      gPBIRStatus;
PBMutateResidue      gPBMRStatus;

PCUpdateConformation   gPCUCStatus;      /* protein conform structures */

HighlightedResidueList   gHighlightedResidues;   /*highlighted residues record*/

short gRemoveTarget;   /* remove molecules target */

/* scan utility record */
ScanRecord gScanStatus;

/* descriptions of scan menu window */
DlogMenuSize ScanMenuSize;
DlogMenuWindow ScanMenuWindow;

/* Display flag pointers & Residue Colour Pointers */
DisplayFlagsCommonPtr gDisplayFlags;
ResidueColoursCommonPtr gResidueColours;

/* Insertion codes pointers */
InsertionCodesCommonPtr gInsertionCodes;

/* Ribbon parameters pointers */
RibbonCommonPtr gRibbon;

/* CGF COMMON block pointer */
ProteinCGFCommonPtr gCGFGroups;

/* Restraints variables */
RestListPtr gRestraintsPtr ;      /* pointer to space for list */
long        gNumRestraints ;      /* number of restraints */
long        gMaxRestraints ;      /* max number of restraints */

/* Tests for whether Asp, Cobra and Constrictor exist */
int         gAspExist ;
int         gCobExist ;
int         gConExist ;
/* OMUPD KYC 12/12/91 add Global for Rattler */
int         gRatExist;

/* Vector list parameters */
VectorListLookupPtr gVectorList;

/* The bond draw list variables */
AtomPtrType   gBondIndex [ MAXATM + MAXDUM + 1 ] ;
BondDraw    * gBondDrawList ;
AtomPtrType * gRibbonIndices ;
short         gNumCtrlRibbon ;
AtomPtrType * gCAOnlyDisplayList ;
AtomPtrType   gNumCAOnlyDisp ;
AtomPtrType   gNumCAOnlyStore ;

/* pointers to COSMIC common blocks */
/* OMUPD KJW 17/05/93 Removed Cosmic globals (CR 1798) */
#if FALSE
CosmicOddsCommon Fodds;
CosmicFlagsCommon Fengsav;
#endif

#endif
