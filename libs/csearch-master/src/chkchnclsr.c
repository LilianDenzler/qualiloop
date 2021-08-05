#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Checks to make sure that a chain closure does not encompass a terminus. */
 
void chkchnclsr(
int startres,
int lastres,
logical *okp
)
{
   int ires;
   int iseg1,iseg2,fstres,lstres;
 
   iseg1 = getseg(&startres,pstruct.segndx,&values.nsegs);
   iseg2 = getseg(&lastres,pstruct.segndx,&values.nsegs);
   fstres = pstruct.segndx[iseg1-1][0] + 1;
   lstres = pstruct.segndx[iseg1][0];
   if (iseg1 != iseg2 || startres <= fstres + 1 || lastres == lstres)
   {
      fprintf(out,"Error in CHECK_CHAIN_CLOSURE -- The chain closure degree \
of freedom\n");
      fprintf(out,"may not span the N or C terminus.\n");
      *okp = false;
   }
   for (ires=startres; ires <= lastres; ires++)
   {
      if (strncmp(pstruct.resnme[ires-1],"PCA ",4) == 0)
         fprintf(out,"Warning from CHECK_CHAIN_CLOSURE -- PCA not handled \
correctly.\n");
   }
}
 
