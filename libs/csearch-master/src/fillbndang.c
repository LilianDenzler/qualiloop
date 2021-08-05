#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Inserts the bond and bond angle in an atom data structure from the
   antecedents specified there as well as the parameters. */
 
void fillbndang(
struct atom *atp
)
{
   if (atp->atomno && atp->ante[2])
   {
      atp->bond = getparbond2(atp->atomno,atp->ante[2]);
      if (atp->ante[1])
         atp->angle = gtprangl2(atp->atomno,atp->ante[2],atp->ante[1]);
   }
}
 
