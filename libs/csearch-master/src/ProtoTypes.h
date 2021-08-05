/******************************************************************************
 *      Name: ProtoTypes.h                                                    *
 *  Function: Declarations of all function.protypes                           *
 * Copyright: (C) Oxford Molecular Limited                                    *
 *----------------------------------------------------------------------------*
 *    Author: G C Calvert Tessella Support Services plc                       *
 *      Date: 28/09/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 * 20-11-90 JPH     Include lots more header files in here, rather than       *
 *                  including them in files already inclded here.             *
 * 05/12/90 GCC     Removed inclusion of *.pro prototype files                *
 * 07/12/90 GCC     Added HP definitions for full ANSII compilation           *
 * 10/12/90 GCC     Added inclusion of Men.h MENU header file                 *
 * 17/01/91 JPH     Added checks for IRIX operating system                    *
 *                   - 3.3 has standard ANSI C , 3.2 does not                 *
 * 21/01/91 GCC     Reordered statement order as the last change caused the   *
 *                  build on the HP to fail.  This is because the definitions *
 *                  used by the HP machine are used within the standard       *
 *                  header files.                                             *
 * 14/05/91 drm     Added protein window fudge factor for DrawProtWin         *
 * 17/07/91 JPH     Added support for IBM RS6000 AIX.                         *
 * 23/07/91 JPH     Seperate out everything into separate logical files.      *
 ******************************************************************************/

#ifndef __PROTOTYPES__
#define __PROTOTYPES__

#include "OMLInclude.h"
#include "PimmsInclude.h"
#include "PGLInclude.h"
#include "EurekaInc.h"

#endif
