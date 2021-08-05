#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Constructs the clump pointed to by clp with free torsion angle, phi. */
 
void constclmp(
struct clump *clp,
float phi
)
{
   float phiat;
   struct atom *atp,**atpp;
 
   clp->cur_phi = phi;
   for (atpp=clp->atoms; (atp = *atpp); atpp++)
   {
      switch (atp->cons_code)
      {
      case st_free:
         phiat = phi;
         break;
      case st_add:
         phiat = phi + atp->offset;
         break;
      case st_fixed:
         phiat = atp->torsion;
         break;
      default:
         fprintf(out,"\nError in constclmp -- Bad cons_code = %d\n",
                atp->cons_code);
         die();
         break;
      }
      constatm(atp,phiat);
   }
}
 
