#include "ProtoTypes.h"
#include "CongenProto.h"
 
void sidechain_f(
struct dof *dofp,
struct sidechain_d *desc
)
{
   switch (desc->sideopt)
   {
   case first:
      /* ACRM added cast */
       sidchnfst(dofp,(struct sidechain_d *)desc,
                      true,desc->residues,(*desc->residues)->clumps);
      break;
   case all:
      /* ACRM added cast */
       sidchnall(dofp,(struct sidechain_d *)desc,
                    desc->residues,(*desc->residues)->clumps);
      break;
   case independent:
      /* ACRM added cast */
      sidchnindi(dofp,(struct sidechain_d *)desc);
      break;
   case combination:
      /* ACRM added cast */
      sidchncomb(dofp,(struct sidechain_d *)desc);
      break;
   case iterative:
   /* ACRM added cast */
      sidchnitr(dofp,(struct sidechain_d *)desc);
      break;
   default:
      fprintf(out,"\nError in SIDECHAIN_F -- Bad sidechain option: %d\n",
             desc->sideopt);
   }
}
 
