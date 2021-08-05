#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* C routine to perform a binary search of an narray looking for numb. */
 
int bin_search(
int numb,
int *narray,
int nlen
)
{
   int istart = 1,
       istop  = nlen,
       nindx;
    do
   {
      nindx = istart + (istop-istart)/2;
 
      if(narray[nindx-1] == numb)
         return(nindx);
      else if(narray[nindx-1] > numb)
         istop = nindx - 1;
      else if(narray[nindx-1] < numb)
         istart = nindx + 1;
   }  while(istop >= istart);
 
   return(0);
}
 
 
