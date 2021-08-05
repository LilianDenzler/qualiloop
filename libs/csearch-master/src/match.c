#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*************************************************************************
match()
Input:      char  *comstring     A character string
            char  *string2       A second string
Output:     int   *nletters      Number of letters matched
Returns:    int                  0 String mismatch
                                 1 First string finished first
                                 2 Second string finished first
 
This routine matches two strings, but stops the comparison as soon
as a space or NULL is found in either string. Th returned value
indicates which string finished first or 0 if the letters before the
space or NULL have a mismatch. The routine calls String ToUpper()
on `comstring' before the comparison.
*************************************************************************/
 
int match(
char comstring[],
char string2[],
int  *nletters
)
{
   int  i;
   char *string1;

   /* V1.2: */
   if(!string2) return(0);

   string1 = (char *)malloc((strlen(comstring) + 2) * sizeof(char));
   struppr(comstring,string1);

   for(i=0;;i++)
   {
      if((!string1[i])||(string1[i]==' '))
      {
         *nletters = i;
         free(string1);
         return(1);
      }
      if((!string2[i])||(string2[i]==' '))
      {
         *nletters = i;
         free(string1);
         return(2);
      }
      if(string1[i] != string2[i])
      {
         *nletters = i;
         free(string1);
         return(0);
      }
   }
   return(0);
}
 
