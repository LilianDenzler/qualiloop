/******************************************************************************
 *      Name: FORTRANTypes                                                    *
 *  Function: Header file to declare types and constants used to communicate  *
 *            with the Nemesis FORTRAN routines.                              *
 *            Each of the structures declared below is an exact analogue of   *
 *            COMMON blocks used inthe FORTRAN part of Pimms.                 *
 *                                                                            *
 *            N.B. The use of the terms, segment and segment pointer, in the  *
 *            context of the structures defined in this header file does not  *
 *            refer to the usual meanings of these terms in C and  Macintosh  *
 *            programming.  In the current context a segment means a segment  *
 *            of a molecule and a segment pointer is an integer which serves  *
 *            as an index into the Atom List.  These terms should not confuse *
 *            since the context will make clear which set of terms is in use. *
 *                                                                            *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Tessella Support Services Ltd.                      *
 *      Date: 24/01/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 * 12/07/90 GCC     Moved NFragmentAtoms to end of structure as this had been *
 *                  causing misalignement of the structure members in the     *
 *                  fortran common block                                      *
 * 18/12/90 GCC     Increased atom limit to 10000                             *
 * 17/01/91 DRM     Updated Fragment limit                                    *
 * 25/01/91 DW      Added display flags structure definitions                 *
 * 12/02/91 DW      Added charge validity definitions                         *
 * 18/02/91 DW      Added CGF COMMON block definitions                        *
 * 26/03/91 GCC     Added new file formats.                                   *
 * 17/07/91 JPH     Changed MAXSEG to MAXSEGMENTS as MAXSEG already defined   *
 *                  in sys/m_param.h on IBM RS6000.                           *
 * 16/10/91 drm     Added definition for invalid charges                      *
 * 01/11/91 DW      Reduce MAXMOL from 999 to 99                              *
 * 04/06/92 KYC     Add set of flags to indicate molecule has been typed      *
 * 17/12/92 KJW     Make sure that PimmsTypes.h is '#include'd                *
 * 26/07/93 KYC     Increase MAXATM to 15000                                  *
 ******************************************************************************/
 
#if !defined(__FORTRANTYPES__)
#define __FORTRANTYPES__

#if !defined(__PIMMSTYPES__)
#include "PimmsTypes.h"
#endif

#define MAXATM     15000   /* maximum number of atoms */
#define MAXMOL        99   /* maximum number of molecules */
#define MAXCON         8   /* maximum number of bonds/atom */
#define MAXSEGMENTS   10   /* maximum number of segments */
#define MAXDUM        10   /* maximum number of dummy atoms */
#define MAXNDC        50   /* maximum components per dummy atom */
#define MAXPLN         5   /* maximum number of planes */
#define MAXNPC         3   /* maximum number of components per plane */
#define MAXNFA       250   /* maximum number of new fragment atoms */

#define NULLSEGMENT   -1   /* -1 indicates a null molecule segment */
 
#define TITLELENGTH   60   /* length of molecule title strings */
#define LABELLENGTH    4   /* length of atomic label strings */
#define RESIDUELENGTH  3   /* length of residue name strings */
#define RESIDUENUMBER  5   /* length of residue number strings */
#define SYMBOLLENGTH   2   /* length of atomic symbol strings */
#define INSERTLENGTH   1   /* Length of insertion codes */
 
/* codes for structure file types */
#define CSSRFile        0
#define PDBFile         1
#define LISTFile        2
#define MOLFile         3
#define SDFile          4
#define COSMICFile      5
#define EUREKAFile      6
#define FDATFile        7

/* Codes for atom label types */
#define ATOM_LABEL     0   /* atom-label label type */
#define ATOM_CHARGE    1   /* partial atomic charge label type */
#define RESIDUE_NAME   2   /* Label by residue name */
#define ATOM_TYPE      3   /* Label by atom type */
 
/* Molecule colour codes */
#define   SPLIT_BOND    0   /* value to flag split bond colouring */
#define   MONO_COLOUR   1   /* value to indicate mono-coloured molecules */

#define   REDINDEX      0   /* Index to red component of colour arrays */
#define   BLUEINDEX     1   /* Index to blue component of colour arrays */
#define   GREENINDEX    2   /* Index to green component of colour arrays */


/* Charge validity indicator */
/* drm 16/10/91  - change valid to 1 and add invalid defn */
#define CHARGES_VALID    1
#define CHARGES_INVALID  0

/* Ribbon parameters */
#define RIBBON_ON        1
#define RIBBON_OFF       0
#define RIBBON_SINGLE    0
#define RIBBON_PROTEIN   1
#define RIBBON_DEFSTRAND 5

#define   ITRUE        -1   /* value of FORTRAN .TRUE. */
#define   IFALSE        0   /* value of FORTRAN .FALSE. */
 
#define ATMCRD       "/ATMCRD/"   /* name of Atom Coords COMMON block */
#define ATMFLG       "/ATMFLG/"   /* name of Atom Flags COMMON block */
#define ATMPRP       "/ATMPRP/"   /* name of Atom Props COMMON block */
#define ATMCON       "/ATMCON/"   /* name of Atom Conects COMMON block */
#define ATMLBL       "/ATMLBL/"   /* name of Atom Label COMMON block */
#define ATMRNO       "/ATMRNO/"   /* name of residue number COMMON block*/
#define BNDDRW       "/BNDDRW/"   /* name of Bond DrawList COMMON block */
#define MOLPTR       "/MOLPTR/"   /* name of Molecule list COMMON block */
#define MOLFLG       "/MOLFLG/"   /* name of Mol flags COMMON block */
#define MOLCOM       "/MOLCOM/"   /* name of Mol centres COMMON block */
#define MOLTTL       "/MOLTTL/"   /* name of Mol titles COMMON block */
#define DUMMY        "/DUMMY/"    /* name of dum at defns COMMON block */
#define PLANES       "/PLANES/"   /* name of plane defs COMMON block */
#define ISOLST       "/ISOLST/"   /* name of isolated atms COMMON block */
#define FRGNEW       "/FRGNEW/"   /* name of New Fragment COMMON block */
 
/* the AtomCoordCommon type is a structure which       */
/* matches the structure of the common block /ATMCRD/  */
typedef struct {
   AtomCoordType x[MAXATM];      /* atomic x - coordinates */
   AtomCoordType y[MAXATM];      /* atomic y - coordinates */
   AtomCoordType z[MAXATM];      /* atomic z - coordinates */
               } AtomCoordCommon, *AtomCoordCommonPtr;
 
/* the AtomPropCommon type is a structure which        */
/* matches the structure of the common block /ATMPRP/  */
typedef struct {
   AtomChargeType charge[MAXATM];            /* atomic partial charges*/
   AtomPtrType    atomicnumber[MAXATM];      /* atomic numbers */
   AtomPtrType    Catomtype[MAXATM];          /* COSMIC Atom type */
               } AtomPropCommon, *AtomPropCommonPtr;
 
/* the AtomFlagsCommon type is a structure which  */
/* matches the structure of the common block /ATMFLG/  */
typedef struct {
   unsigned short   AtomColours[3][MAXATM];      /* atom RGB values */
               } AtomFlagsCommon, *AtomFlagsCommonPtr;
 
/* the AtomConnectCommon type is a structure which        */
/* matches the structure of the common block /ATMCON/  */
typedef struct {
   AtomPtrType connectarray[MAXCON][MAXATM]; /* inter-atomic connectivities */
               } AtomConnectCommon, *AtomConnectCommonPtr;
 
/* the AtomLabelsCommon type is a structure which      */
/* matches the structure of the common block /ATMLBL/  */
typedef struct {
   char   AtomLabel[MAXATM][LABELLENGTH];       /* atomic labels */
   char   AtomResidue[MAXATM][RESIDUELENGTH];   /* atom residue labels */
   short  ResidueNumber[MAXATM];                /* atom residue numbers */
      } AtomLabelsCommon, *AtomLabelsCommonPtr;
 
/* the ResidueNumbersCommon type is a structure which      */
/* matches the structure of the common block /ATMRNO/  */
typedef struct {
   short  IntResidueNumber[MAXATM];   /* internal atomic residue numbers */
      } ResidueNumbersCommon, *ResidueNumbersCommonPtr;
 
/* the MoleculeListCommon type is a structure which   */
/* matches the structure of the common block /MOLPTR/ */
typedef struct
   { short  segmentptr[MAXSEGMENTS][MAXMOL]; /* segment pointers */
     short  segmentlen[MAXSEGMENTS][MAXMOL]; /* segment lengths */
   } MoleculeListCommon, *MoleculeListCommonPtr;
 
/* the BondDrawCommon type is a structure which        */
/* matches the structure of the common block /BNDDRW/  */
typedef struct {
   short  NMolBonds[MAXMOL];              /* number of bonds per molecule */
   AtomPtrType BondDraw[2][2 * MAXATM];   /* list of bonds to draw */
       } BondDrawCommon, *BondDrawCommonPtr;
 
/* the MoleculeFlagsCommon type is a structure which   */
/* matches the structure of the common block /MOLFLG/  */
typedef struct {
   unsigned short qActive[MAXMOL + 1];      /* Molecule Active flags */
   unsigned short qView[MAXMOL + 1];        /* Molecule Visible flags */
   unsigned short qCharge[MAXMOL];          /* valid charges flags */
   unsigned short qHydrogens[MAXMOL];       /* hydrogens visible flags */
   unsigned short qLabels[MAXMOL];          /* labels visible flags */
   unsigned short LabelType[MAXMOL];        /* labels type indicator */
   unsigned short qMolMode[MAXMOL];         /* molecule colour flags */
   unsigned short MolColours[3][MAXMOL + 1];   /* molecule RGB values */
   unsigned short qAutoTyped[MAXMOL];    /* flag indicates mol has been typed*/
      } MoleculeFlagsCommon, *MoleculeFlagsCommonPtr;
 
/* the MoleculeTitlesCommon type is a structure which */
/* matches the structure of the common block /MOLTTL/ */
typedef struct {
   char  MajorTitle[MAXMOL][TITLELENGTH];   /* major Molecule titles */
   char  MinorTitle[MAXMOL][TITLELENGTH];   /* minor Molecule titles */
      } MoleculeTitlesCommon, *MoleculeTitlesCommonPtr;
 
/* the MoleculeCentreCommon type is a structure which */
/* matches the structure of the common block /MOLCOM/ */
typedef struct {
   float x[MAXMOL];      /* molecule centre x - coordinates */
   float y[MAXMOL];      /* molecule centre y - coordinates */
   float z[MAXMOL];      /* molecule centre z - coordinates */
   float xcrot, ycrot, zcrot;   /* centre of rotations */
      } MoleculeCentreCommon, *MoleculeCentreCommonPtr;
 
/* the DummyCommon type is a structure which */
/* matches the structure of the common block /DUMMY/ */
typedef struct {
   AtomPtrType DComponent[MAXNDC][MAXDUM]; /* component atoms of dummy atom */
   short DSize[MAXDUM];         /* size of each dummy atom */
      } DummyCommon, *DummyCommonPtr;
 
/* the PlaneCommon type is a structure which */
/* matches the structure of the common block /PLANES/ */
typedef struct {
   AtomPtrType PComponent[MAXNPC][MAXPLN];   /* component atoms of plane */
      } PlaneCommon, *PlaneCommonPtr;
      
/* the IsolatedAtomCommon type is a structure which */
/* matches the structure of the common block /ISOLST/ */
typedef struct {
   AtomPtrType NIsolatedAtoms;         /* number of isolated atoms */
   AtomPtrType IsolatedAtom[MAXATM];   /* list of isolated atoms */
      } IsolatedAtomCommon, *IsolatedAtomCommonPtr;
      
/* the NewFragmentCommon type is a structure which */
/* matches the structure of the common block /FRGNEW/ */
typedef struct {
   AtomCoordType   x[MAXNFA];       /* atomic x - coordinates */
   AtomCoordType   y[MAXNFA];       /* atomic y - coordinates */
   AtomCoordType   z[MAXNFA];       /* atomic z -coordinates */
   AtomCoordType   xdupl[MAXNFA];   /* duplicate atomic x - coordinates */
   AtomCoordType   ydupl[MAXNFA];   /* duplicate atomic y - coordinates */
   AtomCoordType   zdupl[MAXNFA];   /* duplicate atomic z - coordinates */
   AtomChargeType  echarge[MAXNFA]; /* atomic partial charges*/
   AtomPtrType     atomicnumber[MAXNFA];              /* atomic numbers */
   char            AtomLabel[MAXNFA][LABELLENGTH];    /* atomic labels */
   char            AtomResidue[MAXNFA][RESIDUELENGTH];/* residue labels */
   short           ResidueNumber[MAXNFA];             /* atom residue numbers */
   AtomPtrType     connectarray[MAXCON][MAXNFA];      /* inter-atomic bonds */
   unsigned short  AtomColours[3][MAXNFA];            /* atom RGB values */
   AtomPtrType     NFragmentAtoms;     /* number of atoms in new fragment */
      } NewFragmentCommon, *NewFragmentCommonPtr;

/* Display flags structure and pointer */
typedef struct {
   short dfmain[MAXATM];  /* Main atom display flags */
   short dfbckp[MAXATM];  /* Backup atom display flags */
} DisplayFlagsCommon, *DisplayFlagsCommonPtr;

/* Residue colour structure and pointer */

#define MAX_PROT_DISP_RESIDUES 22

typedef struct {
   short red[MAX_PROT_DISP_RESIDUES];
   short green[MAX_PROT_DISP_RESIDUES];
   short blue[MAX_PROT_DISP_RESIDUES];
} ResidueColoursCommon, *ResidueColoursCommonPtr;

/* Insertion code structure and pointer */

typedef struct {
   char InsertionCode[MAXATM][1];
} InsertionCodesCommon, *InsertionCodesCommonPtr;

/* Ribbon parameters structure */
typedef struct {
   short RibbonState[MAXMOL];     /* Is ribbon on or off ? */
   short RibbonMode[MAXMOL];      /* Colour by underlying protein or single ? */
   short NumberOfStrands[MAXMOL]; /* Number of strands in the ribbon */
   short RibbonRed[MAXMOL];       /* Ribbon red gun value */
   short RibbonGreen[MAXMOL];     /* Ribbon green gun value */
   short RibbonBlue[MAXMOL];      /* Ribbon blue gun value */
} RibbonCommon, *RibbonCommonPtr;

/* CGF groups and colours */

#define CGFLABELLENGTH 80

typedef struct {
   short NumberOfGroups;
   short red[MAX_PROT_DISP_RESIDUES];
   short green[MAX_PROT_DISP_RESIDUES];
   short blue[MAX_PROT_DISP_RESIDUES];
   char label[MAX_PROT_DISP_RESIDUES][CGFLABELLENGTH];
} ProteinCGFCommon, *ProteinCGFCommonPtr;

/* FORTRAN Vector list common block */

typedef struct {
   int NumVectors;
} VectorListLookup, *VectorListLookupPtr;
 
#endif
