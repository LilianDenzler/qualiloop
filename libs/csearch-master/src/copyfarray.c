#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* This is a new routine to replace copy4 in C */
 
void CopyFArray(
float *a,
float *cop,
int   n
)
{
int i;
   for(i=0; i<n; i++)
      cop[i] = a[i];
}
 
