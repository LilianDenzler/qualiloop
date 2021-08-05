/****************************************************************************** 
 *                                                                            *
 *	    Name: CosmicTypes.h                                                   *
 *  Function: CosmicTypes header file - contains definition of the data types *
 *            used by the Nemesis Energy setup and calculation dialogs and    *
 *            windows.  Nemesis uses the COSMIC force field to perform all    *
 *            its energycalculations.                                         *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Tessella Support Services Ltd.                      *
 *      Date: 12/04/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 22/10/90  KJW    Add UPAtomAction field to the EnergyCalcRecord            *
 * 26/05/92  KJW    Updated for C version of COSMIC  					      *
 ******************************************************************************/

#ifndef __COSMICTYPES__
#define __COSMICTYPES__

#ifndef __UNIX__
#ifndef __TYPES__
#include <Types.h>
#endif
#endif

#ifndef __COSMIC__
#include "Cosmic.h"
#endif

/* The EnergyCalcRecord record contains the current */
/* state of the energy calculation settings */
typedef struct 
	{
	Boolean		Optimize;							/* flag for optimized calculations */
	Boolean		OptimizeSwitch[MAX_ENERGY_TERMS];	/* Energy switches for optimized calculations */
	Boolean		SinglePointSwitch[MAX_ENERGY_TERMS];/* Energy switches for single point calculations */
	short		TargetClass;						/* molecule target class */
	short		selectedMolecule;					/* selected molecule */
	short		maxiter;							/* maximumm number of iterations */
	short		UPAtomAction;						/* action when unparameterized atoms found */
	} EnergyCalcRecord;


/* The EnergyTermSwitch record is a structure consisting  */
/* of flags which indicate whether energy terms should be */
/* switched On or OFF during energy calculations          */
typedef struct
   {
   Boolean  qvdw;       /* van der Waals energy term flag */
   Boolean  qelec;      /* electrostatic(Coulomb) Energy term flag */
   Boolean  qbond;      /* bond energy term flag */
   Boolean  qangle;     /* angle energy term flag */
   Boolean  qtorsion;   /* torsional energy term flag */
   } EnergyTermSwitch, *EnergyTermSwitchPtr;

/* The EnergyTerm record is a structure consisting of values  */
/* which are the individual contributions to the total energy */
/* from the different interactions within a molecule.         */
typedef struct
   {
   float Evdw;       /* van der Waals energy term */
   float Eelec;      /* electrostatic(Coulomb) Energy term */
   float Ebond;      /* bond energy term */
   float Eangle;     /* angle energy term */
   float Etorsion;   /* torsional energy term */
   } EnergyTerm, *EnergyTermPtr;

/* the BondParm type is a structure used to  */
/* describe the characteristics of a bond    */
typedef struct
   {
   float   akt;
   float   alo;
   float   bo;
   } bondParm, *bondParmPtr;

/* the angleParm type is a structure used to  */
/* describe the characteristics of an angle   */
typedef struct
   {
   short   ita[3];
   float   aka;
   float   aao;
   } angleParm, *angleParmPtr;

/* the torsionParm type is a structure used to */
/* describe the characteristics of an torsion */
typedef struct
   {
   short   itt[4];
   float   av;
   float   afd;
   } torsionParm, *torsionParmPtr;

/* the vdWParm type is a structure used to describe */
/* the characteristics of each vdW interaction      */
typedef struct
   {
   short   *n3;
   short   *n4;
   float   *rr;
   float   *ee;
   } vdWParm, *vdWParmPtr;

/* the spinBond type is a structure used to identify */
/* the rotateable bonds in a structure               */
typedef struct
   {
   short	nvdind;
   short	rb[4];
   Boolean	lt[MAXATT];
   } spinBond, *spinBondPtr;

/* the PiAtomList type is a structure which */
/* describes all the pi-atoms a molecule    */
typedef struct
   {
   short   piat[MAXATP];
   short   pindex[MAXATT];
   float   hx[MAXATP];
   short   nel[MAXATP];
   } PiAtomList, *PiAtomListPtr;


/* the PiBondList type is a structure which */
/* describes all the pi-bonds in a molecule */
typedef struct
   {
   short   piob[MAXOBP];
   float   kxy[MAXOBP];
   float   sc1[MAXOBP];
   float   sc2[MAXOBP];
   } PiBondList, *PiBondListPtr;

/* the PiSystemList type is a structure which describes all */
/* the atoms and bonds which comprise a single pi-system    */
typedef struct
   {
   short   pisys[MAXATP];
   short   sysat[MXASYS];
   short   sysob[MXBSYS];
   short   sysidx[MAXATP];
   short   pop[MXASYS];
   } PiSystemList, *PiSystemListPtr;

/* the PiSystemSize type is a structure    */
/* which describes the size of a pi system */
typedef struct
   {
   short   natp;
   short   nobp;
   short   natsys;
   short   nobsys;
   short   npsys;
   } PiSystemSize, *PiSystemSizePtr;

/* the CosmicOddsCommon type is a structure which    */
/* matches the structure of the common block /ODDS/  */
typedef struct
   {
   short   nnobs;
   short   natoms;
   short   nvdws;
   short   nangls;
   short   ntors;
   float   sc;
   } CosmicOddsCommon;

/* the CosmicFlagsCommon type is a structure which     */
/* matches the structure of the common block /ENGSAV/  */
typedef struct
   {
   Boolean    used;
   Boolean    chgedf;
   } CosmicFlagsCommon;

/* the HuckCoulomb type is a structure to hold */
/* the Coulomb integral data read from file    */
typedef struct
   {
   short   NumCoulomb;                  /* number of Coulomb integrals */
   short   Type[MAX_HUCK_COULOMB];      /* atom types */
   short   NumPiElec[MAX_HUCK_COULOMB]; /* number of pi electrons */
   float   Integral[MAX_HUCK_COULOMB];  /* the Coulomb integral */
   } HuckCoulomb;

/* the HuckResonance type is a structure to hold */
/* the Resonance integral data read from file    */
typedef struct
   {
   short   NumResonance;                    /* number of Resonance integrals */
   short   iType[MAX_HUCK_RESONANCE];       /* atom types */
   short   jType[MAX_HUCK_RESONANCE];       /* atom types */
   float   Integral[MAX_HUCK_RESONANCE];    /* the Resonance integral */
   float   MinIntegral[MAX_HUCK_RESONANCE]; /* the Resonance integral (min value) */
   float   MaxIntegral[MAX_HUCK_RESONANCE]; /* the Resonance integral (max value) */
   } HuckResonance;

/* the ModBondRecord type is a structure to hold */
/* the bond modification data read from file     */
typedef struct
   {
   short   ptype;         /* atom type */
   short   qtype;         /* atom type */
   float   SLength;       /* single bond length */
   float   SStrength;     /* single bond force constant */
   float   DLength;       /* double bond length */
   float   DStrength;     /* double bond force constant */
   } ModBondRecord;


#endif
