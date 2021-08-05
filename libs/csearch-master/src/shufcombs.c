#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   This recursive routine permutes all the conformations found for each
*   sidechain.
*/
 
void shufcombs(
struct dof *dofp,
struct sidechain_d *desc,
struct sideres **srpp
)
{
   struct sideres *srp;
   struct clump *clp,**clpp,**c2pp;
   int i;
 
   if (srpp == desc->residues) desc->energy = 0.0;
   if ((srp = *srpp) == null)
      dispatch(dofp->next);
   else if (srp->natom == 0)
      shufcombs(dofp,desc,srpp+1);
   else
   {
      for (i=0; i<srp->nconf; i++)
      {
         cpfbest(srp,i);
         for (clpp = srp->clumps; (clp = *clpp); clpp++)
         {
            if (!chkcntcts(clp->atoms,&srp->maxevdw,srp->sidehits))
            {
               if (dbg.cgen)
               {
                  fprintf(out,"Combination failed for residue %.4s %.4s \
Iteration %d\n",pstruct.resnme[srp->resno-1],pstruct.resid[srp->resno-1],i);
               }
               for (c2pp = srp->clumps; c2pp != clpp; c2pp++)
                  delatmsfgrd((*c2pp)->atoms);
               goto nexti;
            }
         }
         if (dbg.cgen > 1)
         {
            fprintf(out,"Combination proceeding for residue %.4s %.4s \
Iteration %d\n",pstruct.resnme[srp->resno-1],pstruct.resid[srp->resno-1],i);
         }
         desc->energy += srp->beste[i];
         shufcombs(dofp,desc,srpp+1);
         desc->energy -= srp->beste[i];
         for (clpp=srp->clumps; (clp = *clpp); clpp++)
            delatmsfgrd(clp->atoms);
nexti:   ;
      }
   }
}
 
