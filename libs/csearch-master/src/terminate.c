
#include "ProtoTypes.h"

#define CR 13
#define LF 10

void terminate(char *string)
{
   int i=0;

   while(string[i])
   {
      if((string[i]==CR)||(string[i]==LF))
      {
         string[i] = '\0';
         break;
      }
      i++;
   }
}
