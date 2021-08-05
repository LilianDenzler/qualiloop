#include "ProtoTypes.h"
#include "CongenProto.h"
 
#ifndef PI
#define PI ((double)4.0*atan((double)1.0))
#endif

/* OMUPD rkw 02/10/92 added RADTODEG to convert from radians to degrees */

#define RADTODEG ((double)180.0 / PI)

/* Calculates NumAngle angles from the atom numbers in At1Array, 
   At2Array and At3Array. Returns the values in AngArray.
   Normally calculates the angle 1-2-3, if the NoSwap flag is not
   set, and the SwapArray flag for this angle is set then the
   angle 1-3-2 is calculated.
   
   02.08.92 Rewritten      Author:  ACRM
*/
 
void angle(short   *At1Array,
           short   *At2Array,
           short   *At3Array,
           logical *SwapArray,
           logical *NoSwap,
           int     *NumAngle,
           single  *AngArray,
           single  *x,
           single  *y,
           single  *z)
{
   single dummy = 9999.0; /* Null coordinate */
   int    AngCount,
          i, j, k;
 
   if(*NumAngle == 0) return;

   for(AngCount=0; AngCount<*NumAngle; AngCount++)
   {
      i = At1Array[AngCount] - 1;
      j = At2Array[AngCount] - 1;
      k = At3Array[AngCount] - 1;

      if(!(*NoSwap))
      {
         if(SwapArray[AngCount])
         {
            k = j;
            j = At3Array[AngCount];
         }
      }
      
      if(i<0 || j<0 || k<0)
      {
         AngArray[AngCount] = 0.0;
      }
      else
      {
         if(x[i]==dummy || x[j]==dummy || x[k]==dummy)
            AngArray[AngCount] = 0.0;
         else
            AngArray[AngCount] = (float)RADTODEG*atomangle(x[i],y[i],z[i],
                                                           x[j],y[j],z[j],
                                                           x[k],y[k],z[k]);
      }
   }
}
 
