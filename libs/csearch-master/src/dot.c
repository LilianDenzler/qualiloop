#include "ProtoTypes.h"
#include "CongenProto.h"
 
float dot(
float *u,
float *v,
int *dim
)
{
int i;
float sum = 0.0;
 
   for(i=0; i < *dim; i++)
      sum += (u[i] * v[i]);
   return(sum);
}
 
