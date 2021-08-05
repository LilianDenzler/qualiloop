#include "ProtoTypes.h"
#include "CongenProto.h"
 
void cgen (
void
)
{
   initfcgen();
   if(pscgncmnd())
   {
      initfgener();
      top_level_env = (jmp_buf *) alloc(sizeof(jmp_buf));
      if(setjmp(*top_level_env) == 0)
         dispatch(dof_head);
      termcgen();
   }
}
 
