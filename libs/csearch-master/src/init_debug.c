#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Initializes the debugger variables to reasonable value. */
 
void init_debug(
void
)
{
   dbg.alloc = 0;	    /* No initialization of alloc returned space */
   dbg.allhp = 1;	    /* Initialize HEAP to 'HEAP' */
   dbg.allstk = 1;	    /* Initialize STACK to 'STCK' */
   dbg.cgen = 0;
   dbg.clschn = 0;
}
 
