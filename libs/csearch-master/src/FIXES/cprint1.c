/* 12.01.94 Fix potential unterminated string bug    By: ACRM
*/

#include "ProtoTypes.h"
#include "CongenProto.h"
 
void cprint1(
char string[],
int *len
)
{
   char outstr[208];
   int i;
 
   if(*len > 200)
   {
      fprintf(out,"Error in CWRITE: String too long\n");
      die();
   }
   for(i=0; i< *len; i++) outstr[i] = string[i];
   outstr[*len] = '\0';

   /* Added ACRM 12.01.94 */
   outstr[199] = '\0';
   
   fprintf(out,"%s\n",outstr);
}
 
