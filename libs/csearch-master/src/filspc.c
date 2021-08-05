#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Fill unused part of string with spaces, doesn't properly terminate
   a C-string */
 
void filspc(
char st[],
int  *stmax,
int *stlen
)
{
int i;
   if(*stlen == *stmax)
      return;
   for(i = *stlen+1; i<*stmax; i++)
      st[i] = ' ';
   return;
}
 
