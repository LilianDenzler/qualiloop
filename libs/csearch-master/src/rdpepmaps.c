#include "ProtoTypes.h"
#include "CongenProto.h"
 
void rdpepmaps(
void
)
{
   /* ACRM changed to C file pointers */
   if (glymap_fp) rdemap(1,&cg.nglymap,&cg.glymap,
                         cg.glyemax,"Glycine");
   if (promap_fp) rdemap(2,&cg.npromap,&cg.promap,
                         cg.proemax,"Proline");
   if (alamap_fp) rdemap(3,&cg.naamap,&cg.aamap,
                         cg.alaemax,"Alanine");
}
 
