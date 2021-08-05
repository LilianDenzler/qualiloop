#include "ProtoTypes.h"
#include "CongenProto.h"

/***************************************************************************/
/*>float atomangle(float xi,float yi,float zi,float xj,float yj,
                   float zj,float xk,float yk,float zk)
   -------------------------------------------------------------
   Input:   float    xi,yi,zi    Input coordinates
                     xj,yj,zj
                     xk,yk,zk
   Returns: float                The angle between the 3 atoms

   Calculates the angle between three coordinates.

   31.07.92 Original.   Author: ACRM
*/
float atomangle(float xi,
                float yi,
                float zi,
                float xj,
                float yj,
                float zj,
                float xk,
                float yk,
                float zk)
{
   float px, py, pz,
         qx, qy, qz,
         sp, sq, 
         a2, cosa2;

   px = xi - xj;
   py = yi - yj;
   pz = zi - zj;
   sp = (float)sqrt((double)(px * px + py * py + pz * pz));

   qx = xk - xj;
   qy = yk - yj;
   qz = zk - zj;
   sq = (float)sqrt((double)(qx * qx + qy * qy + qz * qz));

   cosa2 = (qx * px + qy * py + qz * pz) / (sp * sq);
   a2    = (float)acos((double)cosa2);

   return(a2);
}

