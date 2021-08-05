#include "ProtoTypes.h"
#include "CongenProto.h"
 
#ifdef VAX
#include descrip
#endif
 
void trace(
void
)
{
#ifdef VAX
#ifdef CTRACE
/* Why won't this work??? */
   $DESCRIPTOR(string1," ");
   $DESCRIPTOR(string2,"A traceback will be generated");
   lib$put_output(&string1);
   lib$put_output(&string2);
   lib$signal((ULONG)2);
#else
   fortrace();
#endif
#else
   fprintf(out,"Trace called...\n");
#endif
}
 
