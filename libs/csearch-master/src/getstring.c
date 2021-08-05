#include "ProtoTypes.h"
#include "CongenProto.h"
 
#define DIC 34 /* double inverted commas */
 
int GetString(
char command[],
char strparam[]
)
{
   int i,j,inv_commas;
 
   inv_commas=0;
   j=0;
   for(i=0;;i++)
   {
      if(command[i]==DIC)
      {
         /* Toggle the inv_commas flag */
         inv_commas = abs(inv_commas-1);
         /* Don't copy anything */
         continue;
      }
 
      /* Break out if we're at the end of a line */
      if((command[i]==LF)
       ||(command[i]==CR)
       ||(command[i]=='\0')) break;
 
      /* Also break out if we've a space and we're not between
         inverted commas */
      if((command[i]==' ') && (!inv_commas)) break;
 
      /* Other wise copy the character */
      strparam[j++] = command[i];
   }
   strparam[j]='\0';
   return(i);
}
 
