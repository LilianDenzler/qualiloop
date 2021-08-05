#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Gets the parameter for the bond length between i and j. */
 
float getparbond2(
int i,
int j
)
{
float xxx;

   xxx = mygetparbond(&i,
                      &j,
                      cg.parm_no,
                      pstruct.atcode,
                      engpar.bndkey,
                      &values.nbpar,
                      pstruct.atmnme,
                      engpar.eqbdis);
   return xxx;
}
 
