#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Prints the range list pointed to by rp. */
 
void prtrnglst(
struct range_list *rp
)
{
   float *lowp,*highp;
   int i;
 
   if (rp->nrange == 0)
      fprintf(out,"The list is empty.\n");
   else
   {
      for (i=0, lowp=rp->low, highp=rp->high;
           i < rp->nrange;
           i++,lowp++,highp++)
      {
         fprintf(out,"  %g to %g\n",*lowp,*highp);
      }
   }
}
 
