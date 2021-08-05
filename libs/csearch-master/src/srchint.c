#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Searches array for val. Used in place of srchwd for integers */
 
int srchint(
int *array,
int numval,
int val
)
{
int i;
 
   if(numval == 0) return(0);
   for(i=0; i<numval; i++)
      if(val == array[i]) return(i+1);
   return(0);
}
 
