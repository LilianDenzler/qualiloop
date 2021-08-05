#include "ProtoTypes.h"
#include "CongenProto.h"
 
float vlen(
float *v,
int   *dim
)
{  int i;
   double sumd = 0.0;
   for(i=0; i < *dim; i++)
      sumd += (v[i] * v[i]);
 
   return((float)sqrt(sumd));
}
 
