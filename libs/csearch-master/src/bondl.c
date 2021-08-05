#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Calculates bond lengths of a set of atom number pairs. Normally
   the distance between atoms stored in At1Array and At2Array will
   be calculated. However, if the SwapArray flag is set for this 
   bond, the distance between atoms in At1Array and AltArray will
   be calculated. The bond lengths are returned in BondLengths

   02.08.92 Rewritten   By:   ACRM.
*/

#ifdef DIST
#undef DIST
#endif
#define DIST(a,b) sqrt((double)((x[a] - x[b]) * (x[a] - x[b]) + \
                                (y[a] - y[b]) * (y[a] - y[b]) + \
                                (z[a] - z[b]) * (z[a] - z[b])))

void bondl(short   *At1Array,
           short   *At2Array,
           short   *AltArray,
           logical *SwapArray,
           short   *NumBonds,
           single  *x,
           single  *y,
           single  *z,
           single  *BondLengths)
{
   int    BondCount,
          i, j;
 
   if(*NumBonds == 0) return;

   for(BondCount = 0; BondCount < *NumBonds; BondCount++)
   {
      i = At1Array[BondCount] - 1;
      if(SwapArray[BondCount])   j = AltArray[BondCount] - 1;
      else                       j = At2Array[BondCount] - 1;

      if(i<0 || j<0)
         BondLengths[BondCount] = 0.0;
      else
         BondLengths[BondCount] = (single)DIST(i,j);
   }
}
