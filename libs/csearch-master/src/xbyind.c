#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Exchanges elements in two arrays. The size of the array elements in bytes
*   is given by asize. The arrays are a and b. The array aind is a list of
*   element indices to be swapped, and nind gives the number of indices.
*   The element indices in aind are presumed to start at one, not zero.
*/
 
void xbyind(
int asize,
char a[],
char b[],
int aind[],
int nind
)
{
int i,j;
char t,*ap,*bp;
 
   for (i=0; i<nind; i++)
   {
      ap = a + (aind[i]-1)*asize;
      bp = b + (aind[i]-1)*asize;
      for (j=0; j<asize; j++)
      {
         t = *ap;
         *ap++ = *bp;
         *bp++ = t;
      }
   }
}
 
