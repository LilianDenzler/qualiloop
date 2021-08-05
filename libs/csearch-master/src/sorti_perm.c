#include "ProtoTypes.h"
#include "CongenProto.h"
 
void sorti_perm(
int *array,
int numb,
int *idx
)
{
int *temp,i;
 
   indexi(numb,array,idx);
   temp = (int *)malloc(numb * sizeof(int));
   for(i=0;i<numb;i++)
      temp[i] = array[idx[i]];
   for(i=0;i<numb;i++)
      array[i] = temp[i];
   free(temp);
}
 
