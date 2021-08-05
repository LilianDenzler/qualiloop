#ifdef __FORT_UNDERSCORE__
#define ascale ascale_
#endif

#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Scales an array. Called from FORTRAN so uses pointers

   02.08.92 Some rewriting.   By: ACRM
*/
 
void ascale(float *InArray,
            float *Scale,
            float *OutArray,
            int   *NItems)
{
   int i;
 
   for(i=0; i < *NItems; i++)
      OutArray[i] = *Scale * InArray[i];
}
 
