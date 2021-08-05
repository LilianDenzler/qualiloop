#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Finds sidechain conformations using the independent method. */
 
logical sidchnindi(
struct dof *dofp,
struct sidechain_d *desc
)
{
   struct sideres **srpp,*srp;
 
   desc->energy = 0.0;
   for (srpp = desc->residues; (srp = *srpp); srpp++)
   {
      srp->beste[0] = 0.0;
      if (srp->natom > 0)
      {
         srp->nconf = 0;
         sidchnbst(desc,1,srp,srp->clumps,false);
         if (srp->nconf == 0)
         {
            if (dbg.cgen)
               fprintf(out,"Sidechain_independent: Blocked %.4s %.4s\n",
                      pstruct.resnme[srp->resno-1],pstruct.resid[srp->resno-1]);
            return false;
         }
      }
      desc->energy += srp->beste[0];
   }
   for (srpp = desc->residues; (srp = *srpp); srpp++)
   {
      if (srp->natom > 0)
      {
         cpfbest(srp,0);
         walk_over_sideres(srp,adatmtgrd(atp->atomno););
      }
   }
   dispatch(dofp->next);
   for (srpp = desc->residues; (srp = *srpp); srpp++)
   {
      if (srp->natom > 0)
      walk_over_sideres(srp,delatmfgrd(atp->atomno););
   }
   return true;
}
 
