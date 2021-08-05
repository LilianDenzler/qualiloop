#include "ProtoTypes.h"
#include "CongenProto.h"
 
void struppr(
char string1[],
char string2[]
)
{
int i;
 
   for(i=0;i<strlen(string1);i++)
      string2[i] = toupper(string1[i]);
   string2[i] = '\0';
}
 
