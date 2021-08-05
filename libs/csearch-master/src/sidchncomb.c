#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Performs the combination sidechain placement method. Here, each
*   sidechain is considered independently, and the best ncomb
*   conformations are saved. Then, all the sidechains are combined
*   discarding any combinations that have bad contacts, and for each of
*   these, we dispatch to the next level.
*/
 
logical sidchncomb(
struct dof *dofp,
struct sidechain_d *desc
)
{
   struct sideres **srpp,*srp;
 
   for (srpp = desc->residues; (srp = *srpp); srpp++)
   {
      srp->nconf = 0;
      if (srp->natom > 0)
      {
         sidchnbst(desc,desc->ncomb,srp,srp->clumps,false);
         if (srp->nconf == 0)
         {
            if (dbg.cgen > 0)
               fprintf(out,"Sidechain_combination: Blocked for \
residue %.4s %.4s\n",pstruct.resnme[srp->resno-1],pstruct.resid[srp->resno-1]);
            return false;
         }
      }
   }
   shufcombs(dofp,desc,desc->residues);
   return true;
}
 
