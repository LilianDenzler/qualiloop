#include "ProtoTypes.h"
#include "CongenProto.h"
 
 
/* Given res number, return segment number
   This routine is called from FORTRAN as well as from C,
   so pointers are used in the call */
 
int getseg(
int *ires,
int nictot[maxseg+1][10],
int *nseg
)
{  int start = 0,
       stop  = *nseg - 1,
       try;
 
   do{
      try = (start + stop)/2;
      if(*ires > nictot[try][0] && *ires <= nictot[try+1][0]) return(try+1);
      if(*ires <= nictot[try][0]) stop = try-1;
      if(*ires > nictot[try+1][0]) start = try+1;
   }  while(start <= stop);
   fprintf(out,"Error in GETSEG. Couldn't find residue %d\n",*ires);
   die();
 
   return(0);
}
 
