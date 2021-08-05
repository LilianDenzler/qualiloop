#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Adds the angles present in the chain closure descriptor pointed to
*   by desc into the cangle array.
*/
void adccangs(
int **cangle,
float cmaxdt[],
int *ianglep,
struct chain_closure_d *desc
)
{
   int ir,iangle,i,j,k,t;
 
   for (ir=1; ir<=3; ir++)
   {
      iangle = *ianglep;
      *ianglep = iangle+1;
      i = desc->pos_ca[ir-1];
      j = desc->pos_c[ir-1];
      k = desc->pos_n[ir];
      if (i > k) swap(i,k,t);
      cangle[iangle][0] = i;
      cangle[iangle][1] = j;
      cangle[iangle][2] = k;
      cmaxdt[iangle] = desc->maxdt;
 
      iangle = *ianglep;
      *ianglep = iangle+1;
      i = desc->pos_c[ir-1];
      j = desc->pos_n[ir];
      k = desc->pos_ca[ir];
      if (i > k) swap(i,k,t);
      cangle[iangle][0] = i;
      cangle[iangle][1] = j;
      cangle[iangle][2] = k;
      cmaxdt[iangle] = desc->maxdt;
 
      iangle = *ianglep;
      *ianglep = iangle+1;
      i = desc->pos_n[ir];
      j = desc->pos_ca[ir];
      k = desc->pos_c[ir];
      if (i > k) swap(i,k,t);
      cangle[iangle][0] = i;
      cangle[iangle][1] = j;
      cangle[iangle][2] = k;
      cmaxdt[iangle] = desc->maxdt;
   }
}
 
