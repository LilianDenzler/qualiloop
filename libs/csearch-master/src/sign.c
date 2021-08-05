#include "ProtoTypes.h"
#include "CongenProto.h"

/* OMUPD JAW 11/08/92 Added ProtoTypes.h as part of CR 318 */
 
/*
*   Determines the sign of i.
*/
 
int sign(
int i
)
{
   if (i < 0)
      return -1;
   else if (i == 0)
      return 0;
   else
      return 1;
}
 
