#include "ProtoTypes.h"
#include "CongenProto.h"
 
void padterm(
char string[],
int size
)
{
int i,dopad = 0;
 
   for(i=0;i<size;i++)
   {
      if(dopad)
        string[i] = ' ';
      else if(!string[i])
      {
         dopad = 1;
         string[i] = ' ';
      }
   }
   string[size] = '\0';
}
 
