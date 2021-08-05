/****************************************************************************** 
 *      Name: EurekaGlobal.h                                                  *
 *  Function: Header file to declare global variables for Eureka.             *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: J.P.Holland, Oxford Molecular Ltd.                              *
 *      Date:  7-Sep-1990                                                     *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 14/01/92  KYC    Removed gEurIntPath, gEurProgPath, gEurIntName, gEurProNam*
 *                  Programs now use EnvPath                                  *
 * 10/09/92  JAW    Added new ABM specific globals for CR 520                 *
 ******************************************************************************/
#ifndef __EUREKA_GLOBAL__
#define __EUREKA_GLOBAL__

#include "EurekaConsts.h"

/* Global keyword arrays */
char   *gEurContKeys[EURCD_MAX];   /* Keywords in control deck */
char   *gEurCordKeys[EURCF_MAX];   /* Keywords in coordinate file */
char   *gEurDataKeys[EURDF_MAX];   /* Keywords in data file */
char   *gEurCalcKeys[EURCT_MAX];   /* Keywords in for calculation type */
char   *gEurParmKeys[EURPF_MAX];   /* Keywords in parameter file */
char   *gEurRstrnKeys[EURRF_MAX];  /* Keywords in restraints file */

   /* Various global variables */
/*char   * gEurProgPath ;            /* Path to Eureka program */
/*char   * gEurProgName ;            /* Name of Eureka program */
/*char   * gEurIntPath ;            /* Path to interface program */
/*char   * gEurIntName ;            /* Name of interface program */
char   * gEurControlDeck ;         /* Name of control deck */
char     gEurParamFile [ EURFE_FILENAME_MAX ] ;   /* Eureka param file */
char     gEurekaPath [ PATH_MAX + 1 ] ;      /* Path to output files */
char     gEurInputBody [ NAME_MAX + 1 ] ;
char     gEurOutputBody [ NAME_MAX + 1 ] ;
char     gEurCoordExt [ NAME_MAX + 1 ] ;
char     gEurRestExt [ NAME_MAX + 1 ] ;
char     gEurRstrtExt [ NAME_MAX + 1 ] ;
char     gEurLogExt [ NAME_MAX + 1 ] ;
char     gEurCordInExt [ NAME_MAX + 1 ] ;
char     gEurCordOuExt [ NAME_MAX + 1 ] ;
char     gEurCDExt [ NAME_MAX + 1 ] ;
char     gEurConnExt [ NAME_MAX + 1 ] ;
short    gEurCalcTarget ;

   /* Version number of various files */
int      gEurContVers ;
int      gEurCordVers ;
int      gEurDataVers ;
int      gEurParmVers ;
int      gEurRstrnVers ;

/* Names of months and weekdays */
   /* Names of days of week */
char   * gWeekdayName [ 7 ] ;
   /* Names of months */
char   * gMonthName [ 12 ] ;

/* Calculation type */
short  gEurCalcType ;

long   gEurPimAtom ;
long   gEurCoordSets ;
fpos_t gEurCoordSetPos ;
int    gEurMachineType ;
int    gEurFormatted ;

/* ABM-specific globals follow - copies of these are made in the
   FORTRAN common block ABMEUR defined in the abmeur.inc header file.
   gEurABMMode is set to PIMMS_TRUE if the program has been invoked from 
   ABM, PIMMS_FALSE otherwise.  The setting of the other variables will be
   ignored unless gEurABMMode is PIMMS_TRUE. See abmeur.inc for more
   detailed documentation of the function of these variables. */

long gEurABMMode;    /* Whether running in ABM Mode or not */
long gEurABMCycle;   /* No. of minimisation cycles for ABM */
long gEurABMKeep;    /* No. of minimised structures to keep */
long gEurABMSolvMod; /* Whether to use solvent-modified potential */
char gEurABMOutEn [ EURFE_FILENAME_MAX ];
                     /* Filename for energy output */
char gEurABMOutSt [ EURFE_FILENAME_MAX ];
                     /* Root filename for structure output */
char gEurABMOutEr [ EURFE_FILENAME_MAX ];
                     /* Filename for error output */
char gEurABMInCon [ EURFE_FILENAME_MAX ];
                     /* Filename for conformations input */
char gEurABMPDB   [ EURFE_FILENAME_MAX ];
                     /* Filename for input PDB file - OMUPD JAW 10/09/92 CR 520 */
char gEurABMConv  [ EURFE_FILENAME_MAX ];
                     /* Filename for input atom number conversion file - OMUPD JAW 10/09/92 CR 520 */


#endif
