/****************************************************************************** 
 *      Name: EurekaGlob2.h                                                   *
 *  Function: Declare global variables for Eureka ( in Eureka not Pimms).     *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: J.P.Holland, Oxford Molecular Ltd.                              *
 *      Date: 28-Sep-1990                                                     *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 ******************************************************************************/
#ifndef __EUREKAGLOB2__
#define __EUREKAGLOB2__

	/* Include definition of 'EURFE_FILENAME_MAX' */
#include "EurekaConsts.h"

char    gEurDataFileName  [ EURFE_FILENAME_MAX ] ; /* Store name of datafile  */
char    gEurParamFileName [ EURFE_FILENAME_MAX ] ; /* Store name of param file*/
char    gEurInputCrdName  [ EURFE_FILENAME_MAX ] ; /* Store name of coords in */
char    gEurOutputCrdName [ EURFE_FILENAME_MAX ] ; /* Store name of coord out */
char    gEurRestraintName [ EURFE_FILENAME_MAX ] ; /* Store name of rest file */
char    gEurLogfileName   [ EURFE_FILENAME_MAX ] ; /* Store name of log file */
char    gEurContDeckName  [ EURFE_FILENAME_MAX ] ; /* Store name of cont deck */
char    gEurRestartName   [ EURFE_FILENAME_MAX ] ; /* Store name of restart f */


FILE  * gEurCoordOut ;  /* Stream output file is opened on */
short   gEurStructOut ; /* Structure files to be output every 'n' cycles */
int     gEurPipeOut ;

#endif
