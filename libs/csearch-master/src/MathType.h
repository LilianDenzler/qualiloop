/************************************************************************/
/* Maths type definitions
*/
#ifndef _MATHTYPE_H
#define _MATHTYPE_H

/* This is for compilers running on machines such as Amigas, Macs and
   older Sun workstations using 680X0 series processors with maths
   coprocessors. This assumes that the symbol _M68881 is defined when
   the compiler is run to use the maths coprocessor and that a file
   called m68881.h is to be included to make full use of the coprocessor
*/
#ifdef _M68881
#include <m68881.h>
#endif

/* Note, that if this is changed to float, all I/O routines using type
   REAL will need %lf's changing to %f's
*/
typedef double REAL;

typedef struct
{  REAL x, y, z;
}  VEC3F;

typedef VEC3F COOR;

#endif
