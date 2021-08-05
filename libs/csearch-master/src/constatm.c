#include "ProtoTypes.h"
#include "CongenProto.h"
 
void constatm(
struct atom *atp,
float torsion
)
{
   float phi;
 
   if (atp != null)
   {
      phi = torsion;
      cartx2(coords.xcart,coords.ycart,coords.zcart,
             &atp->ante[0],&atp->ante[1],&atp->ante[2],
             &atp->atomno,&atp->bond,&atp->angle,&phi);
   }
}
 
