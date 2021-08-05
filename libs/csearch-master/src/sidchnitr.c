#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Implements the iterative sidechain placement algorithm. */
 
void sidchnitr(
struct dof *dofp,
struct sidechain_d *desc
)
{
   struct sideres **srpp,*srp,**s2pp,*s2p;
   logical done,none_changed;
   int iter;
   float dele;
 
/* The following case must be handled specially, otherwise, the loop
   below will be infinite. */
 
   if (*dofp->atoms == null)
      dispatch(dofp->next);
   else
   {
      if ( sidchnfst(dofp,desc,false,desc->residues,
         (*desc->residues)->clumps) != 1)
      {
         if(dbg.cgen) fprintf(out,"Sidechain_iterative: Sidechain search failed.\n");
         return;
      }
      for (srpp=desc->residues; (srp = *srpp); srpp++)
      {
         if (srp->natom == 0)
         {
            srp->changed = false;
            srp->beste[1] = 0.0;
         }
         else
         {
            srp->changed = true;
            srp->beste[1] = largnum;
         }
         walk_over_sideres(srp,adatmtgrd(atp->atomno););
         cpintobst(srp,1);
         walk_clumps(srp,clp->best_phi[1] = clp->cur_phi;);
      }
      srpp = desc->residues;
      done = false;
      iter = 0;
      while (!done)
      {
         while ((srp = *srpp)->natom == 0)
            if (*++srpp == null) srpp = desc->residues;
         walk_over_sideres(srp,delatmfgrd(atp->atomno););
         srp->nconf = 0;
         sidchnbst(desc,1,srp,srp->clumps,true);
         iter++;
         if (srp->nconf == 0)
         {
            walk_clumps(srp,constclmp(clp,clp->best_phi[1]);
                        clp->best_phi[0] = clp->best_phi[1]; );
            srp->beste[0] = srp->beste[1];
            srp->nconf = 1;
            cpintobst(srp,1);
            fprintf(out,"Warning from sidchnitr -- \
Sidechain best failed\n");
         }
         if (dbg.cgen)
            fprintf(out,"New sidechain energy for %.4s %.4s is %g. \
Old E is %g.\n",pstruct.resnme[srp->resno-1],pstruct.resid[srp->resno-1],
                       srp->beste[0],srp->beste[1]);
         dele = fabs(srp->beste[0] - srp->beste[1]);
 
         if(dele <= 1.0e-4 ||
            dele / maxf2(fabs(srp->beste[0]),
                         fabs(srp->beste[1])) < 1.0e-4)
         {
            srp->changed = false;
         }
         else
         {
            srp->beste[1] = srp->beste[0];
            cpftbest(srp,0,1);
            srp->changed = true;
            walk_clumps(srp,clp->best_phi[1] = clp->best_phi[0]; );
         }
         cpfbest(srp,1);
         walk_over_sideres(srp,adatmtgrd(atp->atomno););
         none_changed = true;
         for (s2pp = desc->residues; (s2p = *s2pp); s2pp++)
         {
            if (s2p->changed)
            {
               none_changed = false;
               break;
            }
         }
         if (dbg.cgen && none_changed)
            fprintf(out,"Sidechain placed after %d iterations.\n",iter);
         done = none_changed || iter > desc->maxsideiter*desc->nresidues;
         if (*++srpp == null) srpp = desc->residues;
      }
      if (!none_changed)
         fprintf(out,"Sidechain iteration exceeded for conformation %d.\n",
                confnum+1);
      desc->energy = 0.0;
      for (srpp=desc->residues; (srp = *srpp); srpp++)
      {
         desc->energy += srp->beste[1];
      }
      dispatch(dofp->next);
      for (srpp=desc->residues; (srp = *srpp); srpp++)
      {
         walk_over_sideres(srp,delatmfgrd(atp->atomno););
      }
   }
}
 
