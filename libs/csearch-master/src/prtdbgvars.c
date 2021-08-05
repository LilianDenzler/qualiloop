#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Prints the current values for the debugging variables.
*/
 
void prtdbgvars(
void
)
{
   fprintf(out,"\nDebugging variables:\n");
   fprintf(out,"ALLOC = %d  ALLHP = %d  ALLSTK = %d\n",
          dbg.alloc,dbg.allhp,dbg.allstk);
   fprintf(out,"CGEN = %d  CLSCHN = %d\n",dbg.cgen,dbg.clschn);
}
 
