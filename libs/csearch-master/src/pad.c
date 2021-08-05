#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Pad a string with spaces */
 
void pad(
char string[],
int num
)
{
   int i,
       do_rest = 0;
 
   for(i=0;i<num;i++)
   {
      if(do_rest)
      {
         string[i] = ' ';
      }
      else
      {
         if(!string[i])
         {
            do_rest = 1;
            string[i] = ' ';
         }
      }
   }
}
 
