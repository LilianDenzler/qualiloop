#include "ProtoTypes.h"
#include "CongenProto.h"
 
void scnemapln(
char *line,
float *omega,
float *phi,
float *psi,
float *e1,
float *e2
)
{
   sscanf(line,"%10f%10f%10f%11f%15f",omega,phi,psi,e1,e2);
}
 
