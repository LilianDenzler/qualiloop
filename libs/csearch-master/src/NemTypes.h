/******************************************************************************
 *                                                                            *
 *      Name: NemTypes.h                                                      *
 *  Function: Header file to define types and constants for general use in    *
 *            Nemesis. Some of the types defined are intended to make Nemesis *
 *            platform independent.                                           *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Oxford Molecular Limited.                           *
 *      Date: 22/05/92                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 17/11/92  KJW    Made portable across Unix, Macintosh & Windows            *
 * 25/08/93  KJW    Updated for Nemesis v2.0 on the PC                        *
 ******************************************************************************/

#ifndef __NEM_TYPES__
#define __NEM_TYPES__

#ifdef __UNIX__
   #include "TrueFalse.h"
   #include "PimmsTypes.h"
   #ifndef NULL
      #define NULL   0
   #endif
   #define nil      NULL
   #ifndef _XtIntrinsic_h
      #ifndef __GL_GL_H__
         typedef unsigned char   Boolean;   /* define a Boolean type for Unix */
      #endif
   #endif
   typedef void   *H_Handle;  /* define a Handle type for Unix */
#endif

#ifdef __MACINTOSH__
   #include <Types.h>
   #define NULL_CHAR      '\0'
   #define NULL_STR       ""
   enum {FALSE = false, TRUE = true};
   typedef Handle    H_Handle; /* make H_Handle equivalent to the Mac Handle type */
#endif

#ifdef __MS_WINDOWS__
   #ifndef _INC_WINDOWS_
      #include <Windows.h>
   #endif
   typedef BOOL      Boolean;  /* make Boolean equivalent to the Windows BOOL type */
   typedef HANDLE    H_Handle; /* make H_Handle equivalent to the Windows Handle type */
   typedef short     OSErr;    /* define equivalent to Mac OSErr type */
#endif

#define ON       1     /* constants to indicate switch states */
#define OFF      0

#endif  /* __NEM_TYPES__ */
