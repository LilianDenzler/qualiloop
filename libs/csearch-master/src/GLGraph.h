/******************************************************************************
 *      Name: GLGraph.h                                                       *
 *  Function: Header file primarily for use by the PGL library when           *
 *            used on the machines using the GL library.                      *
 * Copyright: (C) Oxford Molecular Limited                                    *
 *----------------------------------------------------------------------------*
 *    Author: J.P.Holland, Oxford Molecular Ltd.                              *
 *      Date: 26-Jul-1991                                                     *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 28/12/91 drm     Added the font manager interface for the Silicon Graphics *
 ******************************************************************************/

#ifndef __GLGRAPH__
#define __GLGRAPH__

#if defined (__SILICONGRAPHICS__)
#include <gl.h>
#include <gl/get.h>
#include <device.h>
/* drm 28/12/91 */
#include <fmclient.h>
#else						/* NOT __SILICONGRAPHICS__ */

#if defined (__IBM__)
#include <gl.h>
#include <device.h>
#else						/* NOT __IBM__ */

/* Define tgl types used in Pimms for other graphics libraries */
typedef unsigned short Colorindex ;
typedef unsigned short Linestyle  ;
typedef unsigned short Device     ;
typedef float          Coord      ;
typedef float          Matrix [ 4 ] [ 4 ] ;
typedef short Screencoord ;

#endif						/* END __IBM__ */
#endif						/* END __SILICONGRAPHICS__ */

#endif	/* __GLGRAPH__ */


