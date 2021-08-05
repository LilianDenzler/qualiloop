#include "ProtoTypes.h"
#include "CongenProto.h"
 
void clear_atom(
struct atom *p
)
{
   int j;
 
   for (j=0; j<=2; j++) p->ante[j] = 0;
   p->bond = 0.0;
   p->angle = 0.0;
   p->torsion = 0.0;
   p->torsion_period = 0.0;
   p->offset = 0.0;
   p->cons_code = 0;
}
 
