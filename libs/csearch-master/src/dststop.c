#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Branches (using longjmp) out of all dispatch calls. */
 
void dststop(
void
)
{
   fprintf(out,"\nCGEN execution terminating prior to completion of search\n");
   longjmp(*top_level_env,1);
}
 
