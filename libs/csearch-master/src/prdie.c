#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Prints a message then calls die() */
 
void prdie(
char string[]
)
{
   printf("%s",string);
   die();
}
 
