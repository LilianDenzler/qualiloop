#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* This does a rdtitl. It echoes the title lines to the output list file. */
 
void rdprtitle(
FILE *fp
)
{
   char buffer[201];
 
   for(;;)
   {
      if(!fgets(buffer,200,fp))
         break;
      terminate(buffer);
      if(buffer[0] == '*')
      {
          fprintf(out,"\n");
          break;
      }
 
      fprintf(stdout,"%s\n",buffer);
   }
}
 
