#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Evaluates the current coordinates by computing the rms deviations
*   from the reference coordinates.
*/
 
double sidchnrms(
struct sideres *srp
)
{
   double sum2;
   int *ip,i;
   float *xp,*yp,*zp,dx,dy,dz;
 
   xp = srp->refx;
   yp = srp->refy;
   zp = srp->refz;
   sum2 = 0.0;
   for (ip = srp->atomp; (i = *ip); ip++)
   {
      dx = coords.xcart[--i] - *xp++;
      dy = coords.ycart[i] - *yp++;
      dz = coords.zcart[i] - *zp++;
      sum2 += dx*dx + dy*dy + dz*dz;
   }
   return srp->natom == 0 ? 0.0 : sqrt(sum2/srp->natom);
}
 
