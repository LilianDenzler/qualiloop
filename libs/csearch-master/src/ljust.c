#include "ProtoTypes.h"
#include "CongenProto.h"
 
void ljust(
char string[]
)
{
   int len, j, k;
   char *temp;
 
   len = strlen(string);
   temp = (char *)malloc((len+1) * sizeof(char));
 
   j=0; k=0;
   while(string[j])
   {
      if(string[j] != ' ') temp[k++] = string[j];
      j++;
   }
   temp[k] = '\0';
   strcpy(string,temp);
   free(temp);
 }
 
