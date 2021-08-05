#include "ProtoTypes.h"
#include "CongenProto.h"
 
void ljustpad(
char string[]
)
{
   ljust(string);
   padterm(string,4);
}
 
