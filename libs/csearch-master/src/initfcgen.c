#include "ProtoTypes.h"
#include "CongenProto.h"
 
 
void initfcgen (
void
)
{
   dof_head = null;
   dof_tail = null;
   setup_parmno(cg.parm_no);
   confnum = 0;
   leafnum = 0;
   cg.maxdt_def = 5.0;
}
 
