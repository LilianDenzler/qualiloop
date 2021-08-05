#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* This routine strips leading spaces from a string. */
 
char *killspcs(
char string[]
)
{
   while(*string == ' ') string++;
   return(string);
}
 
