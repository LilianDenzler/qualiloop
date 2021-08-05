#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Sets up imove for the evaluation degree of freedom. */
 
void setup_imove(
struct evaluate_d *desc
)
{
int i;
short i1;
 
   desc->imove = (short int *) alloc(values.natoms * sizeof(short int));
   /* ACRM Should be no problems here; added cast */
   fill2(desc->imove,&values.natoms,(i1=(short)1,&i1));
   for (i=0; i<desc->nconsp; i++)
     desc->imove[desc->consp[i]-1] = 0;
}
 
