#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Copies the current coordinates into the best coordinate array in
*   in srp. ind gives the index in the array where the index starts
*   at 1. */
 
void cpintobst(
struct sideres *srp,
int ind
)
{
   int *ip,i;
   float *xp,*yp,*zp;
 
   xp = srp->bestxpp[ind];
   yp = srp->bestypp[ind];
   zp = srp->bestzpp[ind];
   for (ip = srp->atomp; (i = *ip); ip++)
   {
      *xp++ = coords.xcart[--i];
      *yp++ = coords.ycart[i];
      *zp++ = coords.zcart[i];
   }
}
 
