#include "ProtoTypes.h"
#include "CongenProto.h"
 
char *alloc(
int n
)
{
   char *p;
   int i;

   if(n == 0)
   {
      p = NULL;
      return p;
   }
 
   p = (char *)malloc(n);
   if(!p)
   {
      fprintf(out," Error in ALLOC -- Insufficient memory\n");
      die();
   }
 
   if (dbg.alloc > 0)
   {
      for (i=0; i<n; i++) p[i] = -1;
   }
   return p;
}
 
