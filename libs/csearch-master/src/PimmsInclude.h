/******************************************************************************
 *      Name: PimmsInclude.h                                                  *
 *  Function: All main Pimms include files.                                   *
 *----------------------------------------------------------------------------*
 *    Author: J.P.Holland, Oxford Molecular Ltd.                              *
 *      Date: 26-Jul-1991                                                     *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 ******************************************************************************/

#ifndef __PIMMSINCLUDE__
#define __PIMMSINCLUDE__

#include "FORTRANCalls.h"
#include "PimmsTypes.h"
#include "Colours.h"
#include "Graph3D.h"
#include "FORTRANTypes.h"

#include "ActiveMols.h"
#include "DispSetup.h"
#include "DispSetupTyp.h"
#include "FileIO.h"
#include "Restraints.h"
#include "RingDef.h"
#include "UtilityProto.h"

#include "ConformTyp.h"		/* Must be AFTER UtilityProto.h */

#ifdef __HEWLETTPACKARD__
#define	PROT_WIN_LENGTH_FUDGE	0
#define	PROT_WIN_HEIGHT_FUDGE	0
#else
#define	PROT_WIN_LENGTH_FUDGE	-1
#define	PROT_WIN_HEIGHT_FUDGE	-1
#endif

#endif
