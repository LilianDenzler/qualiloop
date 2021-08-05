#include "ProtoTypes.h"
#include "CongenProto.h"
 
void fill4(
long *a,
int *n,
long *value
)
{
int i;
   for(i=0; i< *n; i++)
      a[i] = *value;
}
 
