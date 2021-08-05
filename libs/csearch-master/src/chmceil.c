#include "ProtoTypes.h"
#include "CongenProto.h"
 
int chmceil(
float *r
)
{
   int retval;
 
   if(*r < 0.0)
      retval = (int)(*r);
   else if(*r == 0.0)
      retval = 0;
   else
   {
      retval = (int)(*r);
      if(*r - (float)retval != 0.0) retval++;
   }
   return(retval);
}
