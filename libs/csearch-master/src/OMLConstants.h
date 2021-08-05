/*******************************************************************************
**        Name: OMLConstants.h                                                 *
**    Function: Mathematics and other useful constants. Blatantly ripped off   *
**              from math.h in parts                                           *
**   Copyright: (C) Oxford Molecular Limited                                   *
**-----------------------------------------------------------------------------*
**      Author: Matthew Seaman, Oxford Molecular Limited.                      *
**        Date: 23-Jul-1991                                                    *
**-----------------------------------------------------------------------------*
**-----------------------------------------------------------------------------*
** Modification Record                                                         *
** Date     Inits    Comments                                                  *
** dd/mm/yy                                                                    *
*******************************************************************************/

#ifndef __OMLCONSTANTS_H__
#define __OMLCONSTANTS_H__

/*
** Some useful Mathematical constants
*/

#ifdef	M_E
#undef	M_E
#endif
#define	M_E		2.71828182845904523536

#ifdef	M_LOG2E
#undef	M_LOG2E
#endif
#define	M_LOG2E		1.44269504088896340736

#ifdef	M_LOG10E
#undef	M_LOG10E
#endif
#define	M_LOG10E	0.43429448190325182765

#ifdef	M_LN2
#undef	M_LN2
#endif
#define	M_LN2		0.69314718055994530942

#ifdef	M_LN10
#undef	M_LN10
#endif
#define M_LN10		2.30258509299404568402

#ifdef	M_PI
#undef	M_PI
#endif
#define	M_PI		3.14159265358979323846

#ifdef	M_PI_2
#undef	M_PI_2
#endif
#define	M_PI_2		1.57079632679489661923

#ifdef	M_PI_4
#undef	M_PI_4
#endif
#define	M_PI_4		0.78539816339744830962

#ifdef	M_1_PI
#undef	M_1_PI
#endif
#define	M_1_PI		0.31830988618379067154

#ifdef	M_2_PI
#undef	M_2_PI
#endif
#define M_2_PI		0.63661977236758134308

#ifdef	M_2_SQRTPI
#undef	M_2_SQRTPI
#endif
#define	M_2_SQRTPI	1.12837916709551257390

#ifdef	M_SQRT2
#undef	M_SQRT2
#endif
#define	M_SQRT2		1.41421356237309504880

#ifdef	M_SQRT1_2
#undef	M_SQRT1_2
#endif
#define M_SQRT1_2	0.70710678118654752440

/*
** Conversion factors: radians <-> degrees
*/

#ifdef RAD_TO_DEG
#undef RAD_TO_DEG
#endif
#define RAD_TO_DEG	57.2957795130823208768

#ifdef DEG_TO_RAD
#undef DEG_TO_RAD
#endif
#define DEG_TO_RAD	0.01745329251994329576

#ifdef rad120
#undef rad120
#endif
#define rad120 (120.0/180.0*3.141592653589794)
#define dtorad (3.141592653589794/180.0)
#define largnum 1.7e38
#define anum 9999.0
#define largint 2147483647
#define eps1  1.0e-6

/*
** Physical Constants
*/

#endif /* __OMLCONSTANTS_H__ */
