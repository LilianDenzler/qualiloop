#include "ProtoTypes.h"
#include "CongenProto.h"
 
float mincoor(
float *a,
int   n
)
{
   float minval = 9999.0;
   int   i;
 
   for(i=0;i<n;i++)
      if(a[i] != anum) minval = minf2(minval,a[i]);
   if(minval == anum) minval = 0.0;
 
   return(minval);
}
 
