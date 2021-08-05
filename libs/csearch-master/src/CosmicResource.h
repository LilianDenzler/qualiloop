/****************************************************************************** 
 *                                                                            *
 *	    Name: CosmicResource.h                                                *
 *  Function: Cosmic header file - contains definition of constants used by   *
 *            the Resource compiler to build Energy calculation dialogs.      *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Tessella Support Services Ltd.                      *
 *      Date: 28/05/92                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy         													      *
 ******************************************************************************/

#ifndef __COSMIC_RESOURCE__
#define __COSMIC_RESOURCE__

/* Energy Calculation distance dialog item identifiers */
#define dEnergySetupDlg			1501		/* energy setup dialog resource ID */
#define dEnergySetupDitl		1501		/* energy setup dialog item list resource ID */
#define iESSetupParamButton		   3		/* energy setup dialog item list item ID */
#define iESAllActive			   5		/* energy setup dialog item list item ID */
#define iESAllVisible 			   6		/* energy setup dialog item list item ID */
#define iESAllMolecules			   7		/* energy setup dialog item list item ID */
#define iESOneMolecule			   8		/* energy setup dialog item list item ID */
#define iESOneNumber			  10		/* energy setup dialog item list item ID */
#define iESOptimSwitch			  18		/* energy setup dialog item list item ID */
#define iESOptimizationOff		  11		/* energy setup dialog item list item ID */
#define iESOptimizationOn		  12		/* energy setup dialog item list item ID */
#define iESIterations			  14		/* energy setup dialog item list item ID */
#define iESUPAtomAbortOn		  16		/* energy setup dialog item list item ID */
#define iESUPAtomAbortOff		  17		/* energy setup dialog item list item ID */
#define iESHighlight			  18		/* energy setup dialog item list item ID */

/* Energy Calculation distance dialog item identifiers */
#define dEnergyParametersDlg	1502		/* energy setup dialog resource ID */
#define dEnergyParametersDitl	1502		/* energy setup dialog item list resource ID */
#define iEPSPSwitch				   5		/* energy setup dialog item list item ID */
#define iEPOptimSwitch			  11		/* energy setup dialog item list item ID */
#define iEPHighlight			  13		/* energy setup dialog item list item ID */
#define iEPvanderWaals			   0		/* energy setup - van der Waals switch offset */
#define iEPElectrostatics		   1		/* energy setup - electrostatics switch offset */
#define iEPBonds				   2		/* energy setup - bonds switch offset */
#define iEPAngles				   3		/* energy setup - angles switch offset */
#define iEPTorsions				   4		/* energy setup - torsions switch offset */

/* COSMIC Types error dialog item identifiers */
#define dTypesErrorDlg			1503		/* COSMIC Types error dialog resource ID */
#define dTypesErrorDitl			1503		/* COSMIC Types error dialog item list resource ID */
#define iTEHighlight			   6		/* COSMIC Types error dialog item list item ID */

#endif
