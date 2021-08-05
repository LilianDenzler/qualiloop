#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Checks that the range of residues for backbone degree of freedom
   does not span N termini or C termini nor does it hit unplanned for
   residues.*/
 
void chkbbone(
int startres,
int lastres,
logical *okp
)
{
   int ires;
 
   for (ires=startres; ires <= lastres; ires++)
   {
      if(strncmp(pstruct.resnme[ires-1],"NTER",4) == 0 ||
         strncmp(pstruct.resnme[ires-1],"CTER",4) == 0)
      {
         fprintf(out,"Error in CHECK_BACKBONE -- The backbone degree of freedom\n");
         fprintf(out,"may not span the N or C terminus.\n");
         *okp = false;
      }
      if (strncmp(pstruct.resnme[ires-1],"PCA ",4) == 0)
         fprintf(out,"Warning from CHECK_BACKBONE -- PCA not handled \
correctly.\n");
   }
}
 
