#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* A variation of filspc which works with C strings.
   It takes the current string and, if necessary, adds trailing
   spaces to pad to length with a NULL afterwards. */
 
void padspace(
char string[],
int  length
)
{
int i,current_len;
 
   current_len = strlen(string);
   if(current_len < length)
   {
      for(i=current_len; i<length; i++) string[i] = ' ';
      string[length] = '\0';
   }
}
 
