#include "ProtoTypes.h"
#include "CongenProto.h"
 
void fill2(
short *a,
int *n,
short *value
)
{
int i;
 
   for(i=0; i< *n; i++)
      a[i] = *value;
}
 
