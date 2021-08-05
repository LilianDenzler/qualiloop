#include "ProtoTypes.h"
#include "CongenProto.h"
 
void skiptitle(
FILE *fp
)
{
char buffer[201];
 
   for(;;)
   {
      if(!fgets(buffer,200,fp)) break;
      terminate(buffer);
      if(buffer[0] == '*') break;
   }
}
 
