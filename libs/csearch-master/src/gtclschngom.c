#include "ProtoTypes.h"
#include "CongenProto.h"
 
 
/* Get the geometric parameters for use by clschn. */
 
void gtclschngom(
struct chain_closure_d *dp
)
{
   int ir;
 
   for (ir=1; ir<=3; ir++)
   {
      dp->clsbond[ir-1][0] = getclsbond(dp->pos_ca[ir-1],dp->pos_c[ir-1]);
      dp->clsbond[ir-1][1] = getclsbond(dp->pos_c[ir-1],dp->pos_n[ir]);
      dp->clsbond[ir-1][2] = getclsbond(dp->pos_n[ir],dp->pos_ca[ir]);
      dp->clsangle[ir-1][0] = getclsangle(dp->pos_ca[ir-1],dp->pos_c[ir-1],
                                          dp->pos_n[ir]);
      dp->clsangle[ir-1][1] = getclsangle(dp->pos_c[ir-1],dp->pos_n[ir],
                                          dp->pos_ca[ir]);
      dp->clsangle[ir-1][2] = getclsangle(dp->pos_n[ir],dp->pos_ca[ir],
                                          dp->pos_c[ir]);
   }
}
 
