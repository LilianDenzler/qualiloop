#include "ProtoTypes.h"
#include "CongenProto.h"

 
int gcd(
int in,
int jn
)
{
int i,j,k;
 
   if(abs(in) > abs(jn))
   {
      i = abs(in);
      k = abs(jn);
   }
   else
   {
      k = abs(in);
      i = abs(jn);
   }
 
   do {
      j = i;
      i = k;
      k = j%i;
   } while(k!=0);
 
   return(i);
}
 
