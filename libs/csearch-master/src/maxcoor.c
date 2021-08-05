#include "ProtoTypes.h"
#include "CongenProto.h"
 
float maxcoor(
float *a,
int   n
)
{
   float maxval = -9999.0;
   int   i;
 
   for(i=0;i<n;i++)
      if(a[i] != anum) maxval = maxf2(maxval,a[i]);
   if(maxval == -anum) maxval = 0.0;
 
   return(maxval);
}
 
