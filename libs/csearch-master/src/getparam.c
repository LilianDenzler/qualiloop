#include "ProtoTypes.h"
#include "CongenProto.h"
 
int GetParam(
char  command[],
float *value,
int   *nletters
)
{
   char buffer[50];
   int  retval;
 
   if((*nletters = GetString(command,buffer))==0)
      return(0);
 
   retval = sscanf(buffer,"%f",value);
   return(retval);
}
 
