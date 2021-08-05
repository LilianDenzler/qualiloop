#include "ProtoTypes.h"
#include "CongenProto.h"
 
 
/*
*   Add all atoms in the list to the close contact grid.
*/
void adatmstgrd(
struct atom **atoms
)
{
   struct atom *atp;
 
   for (;(atp = *atoms); atoms++) adatmtgrd(atp->atomno);
}
 
