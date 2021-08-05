#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Copies the current coordinates from the best coordinate array in
*   in srp. ind gives the index in the array where the index starts
*   at 1. */
 
void cpfbest(
struct sideres *srp,
int ind
)
{
   int *ip,i;
   float *xp,*yp,*zp;
 
   xp = srp->bestxpp[ind];
   yp = srp->bestypp[ind];
   zp = srp->bestzpp[ind];
   for (ip = srp->atomp; *ip != 0 ; ip++)
   {
      i = ( * ip ) - 1 ;
      coords.xcart[i] = *xp++;
      coords.ycart[i] = *yp++;
      coords.zcart[i] = *zp++;
   }
}
 
