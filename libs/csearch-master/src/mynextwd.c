#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* ACRM: gets a space, tab or null delimited word out of a string */
 
#define string (*instring)
 
void mynextwd(
char *instring[],
char word[]
)
{
   int i;
 
   /* Remove leading spaces */
   while(*string==' ' || *string=='\t') string++;
 
   /* Copy the space or null delimited string into word */
   for(i=0; *string!=' ' && *string!='\t' && *string; i++)
   {
      word[i] = *string;
      string++;
   }
 
   /* Null terminate */
   word[i] = '\0';
 }
 
