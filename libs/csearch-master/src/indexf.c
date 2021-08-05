#include "ProtoTypes.h"
#include "CongenProto.h"
 
/***************************************************************
      Andrew C.R. Martin                            23.06.90
 
      Laboratory of Mathematical Biology            V1.0
      National Institue for Medical Research,
      The Ridgeway,
      Mill Hill,
      London,
      NW7 1AA
****************************************************************
   indexf(n,arrin,indx)
   Input:  n      int      Number of elements in array
           arrin  *float   Array to be indexed
   Output: indx   *int     Index array
 
   This routine uses a heapsort to index a floating point array
   such that arrin[indx[j]] is in ascending order with j.
   It is modified from the FORTRAN version in 'Numerical Recipes'
   Page 233. This version correctly sorts from array element 0
   as opposed to 1 in the FORTRAN version.
****************************************************************/
 
void indexf(
int n,
float *arrin,
int   *indx
)
{
   int   i,j,l,ir,indxt;
   float q;
 
/* OMUPD rkw 04/06/93 Only sort if there are more than 1 items */

   if (n==1)
   {
      indx[0] = 0;
      return;
   }


   for(j=0; j<n; j++) indx[j]=j;
 
   l=n/2+1;
   ir=n;
 
   while(1)
   {
      if(l>1)
      {
         indxt=indx[--l - 1];
         q=arrin[indxt];
      }
      else
      {
         indxt=indx[ir-1];
         q=arrin[indxt];
         indx[ir-1]=indx[0];
         if(--ir == 1)
         {
            indx[0]=indxt;
            return;
         }
      }
      i=l;
      j=l+1;
 
      while(j<=ir)
      {
         if(j<ir)
         {
            if(arrin[indx[j-1]] < arrin[indx[j]]) j++;
         }
         if(q<arrin[indx[j-1]])
         {
            indx[i-1]=indx[j-1];
            i=j;
            j+=j;
         }
         else
         {
            j=ir+1;
         }
      }
      indx[i-1]=indxt;
   }
   return;
}
 
