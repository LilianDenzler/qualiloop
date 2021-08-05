#include "ProtoTypes.h"
#include "CongenProto.h"
 
void dispatch(
struct dof *dofp
)
{
   if (dofp == null)
   {
      if (++leafnum >= cg.maxleaf)
      {
         fprintf(out,"\nMaxleaf has been exceeded.\n");
         dststop();
      }
      return;
   }
   switch(dofp->dof_type)
   {
   case chain_closure_t:
      /* ACRM added cast */
      chnclsr_f(dofp,(struct chain_closure_d *)dofp->desc);
      break;
   case backbone_t:
      /* ACRM added cast */
      backbone_f(dofp,(struct backbone_d *)dofp->desc);
      break;
   case sidechain_t:
      /* ACRM added cast */
      sidechain_f(dofp,(struct sidechain_d *)dofp->desc);
      break;
   case write_coordinates_t:
      /* ACRM added cast */
      wrtcoordsf(dofp,(struct write_coordinates_d *)dofp->desc);
      break;
   case status_t:
      /* ACRM added cast */
      status_f(dofp,(struct status_d *)dofp->desc);
      break;
   case rbest_t:
      /* ACRM added cast */
      rbest_f(dofp,(struct rbest_d *)dofp->desc);
      break;
#ifdef EVALUATE
   case evaluate_t:
      /* ACRM added cast */
      evaluate_f(dofp,(struct evaluate_d *)dofp->desc);
      break;
#endif /* EVALUATE */
   default:
      fprintf(out,"Error in DISPATCH -- Unknown type code %d\n",
                 dofp->dof_type);
      die();
   }
}
 
