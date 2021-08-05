#include "ProtoTypes.h"
#include "CongenProto.h" 
 
/* This routine actually calculates the CTER OT2 coords. */
 
int ccterxyz (
PDB *p,      /* Pointer to PDB item for OT2           */
PDB *ca_p,   /* Pointer to PDB item for anticedant CA */
PDB *c_p,    /* Pointer to PDB for anticedant C       */
PDB *o_p     /* Pointer to PDB for anticedant OT1     */
)
{
   float fac,
         gr = 1.3,
         alpha = 120.0*PI/180.0,
         cosa,  sina,  scalpr,
         x21,   y21,   z21,   d21,
         x23,   y23,   z23,   d23,
         x32,   y32,   z32,
         xp23,  yp23,  zp23,  rp23,
         xh,    yh,    zh,
         xp,    yp,    zp,
         xv,    yv,    zv;
 
   if(!ca_p || !c_p || !o_p)
   {
      printf("ccterxyz(): Error==> Atom pointer undefined.\n");
      return(0);
   }
 
   fac = 0.5 * (float)sqrt((double)3.0);
 
   x21 = c_p->x - ca_p->x;
   y21 = c_p->y - ca_p->y;
   z21 = c_p->z - ca_p->z;
   d21 = (float)sqrt((double)(x21*x21 + y21*y21 + z21*z21));
 
   x23 = c_p->x - o_p->x;
   y23 = c_p->y - o_p->y;
   z23 = c_p->z - o_p->z;
   d23 = (float)sqrt((double)(x23*x23 + y23*y23 + z23*z23));
 
   cosa = (float)cos((double)alpha);
   sina = (float)sin((double)alpha);
 
   x32 = o_p->x - c_p->x;
   y32 = o_p->y - c_p->y;
   z32 = o_p->z - c_p->z;
 
   scalpr = (x21*x32 + y21*y32 + z21*z32)/d21;
   xh = x21 / d21;
   yh = y21 / d21;
   zh = z21 / d21;
 
   xp = scalpr * xh;
   yp = scalpr * yh;
   zp = scalpr * zh;
 
   xp23 = xp + x23;
   yp23 = yp + y23;
   zp23 = zp + z23;
 
   rp23 = (float)sqrt((double)(xp23*xp23 + yp23*yp23 + zp23*zp23));
   xv = xp23/rp23;
   yv = yp23/rp23;
   zv = zp23/rp23;
 
   p->x = c_p->x + gr*(-cosa*xh + sina*xv);
   p->y = c_p->y + gr*(-cosa*yh + sina*yv);
   p->z = c_p->z + gr*(-cosa*zh + sina*zv);
 
   return(0);
}
 
