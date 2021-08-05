#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Gets the parameter for the bond angle between i, j, and k. */
 
float gtprangl2(
int i,
int j,
int k
)
{
   float xxx;

   /* xxx = gtprangl(&i,&j,&k,cg.parm_no,pstruct.atcode,engpar.angkey,
                     &values.napar,&values.natyps,pstruct.atmnme,engpar.eqang);
   */

   dwprangl(&i,&j,&k,
            cg.parm_no,
            pstruct.atcode, engpar.angkey, &values.napar, &values.natyps,
            pstruct.atmnme,
            engpar.eqang, &xxx);

   return xxx;
}
 
