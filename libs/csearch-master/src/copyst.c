#include "ProtoTypes.h"
#include "CongenProto.h"
 
void copyst(
char st1[],
int *st1siz,
int *st1len,
char st2[],
int *st2len
)
{
int lim,i;
 
   if(*st2len == 0)
   {
      *st1len = 0;
      return;
   }
   if(*st1siz <= 0)
   {
      fprintf(out,"Warning from COPYST - Destination string to small\n");
      trace();
   }
   else
   {
      lim = minf2(*st1siz, *st2len);
      for(i=0;i<lim;i++) st1[i] = st2[i];
      *st1len = lim;
      if(*st1siz > *st2len)
      {
         fprintf(out,"Warning from COPYST - Destination string to small\n");
         trace();
      }
   }
}
 
