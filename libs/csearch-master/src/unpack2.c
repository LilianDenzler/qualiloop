#include "ProtoTypes.h"
#include "CongenProto.h"
 
void unpack2(
unsigned long int i4,
unsigned short int *i2,
unsigned short int *j2
)
{
   *i2 = (i4+(1<<15))/(1<<16);
   *j2 = i4 - *i2*(1<<16);
}
 
