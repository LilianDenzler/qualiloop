#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Delete all atoms in the list from the close contact grid. */
 
void delatmsfgrd(
struct atom **atoms
)
{
struct atom *atp;
   for (;(atp = *atoms); atoms++)
      delatmfgrd(atp->atomno);
}
 
