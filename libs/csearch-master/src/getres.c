#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Given atom number, return res number */
 
int getres(
short atom,
short *ibase,
int nres
)
{  int start = 0,
       stop  = nres-1,
       try;
 
   do{
      try = (start + stop)/2;
      if(atom >  ibase[try] && atom <= ibase[try+1]) return(try+1);
      if(atom <= ibase[try]) stop = try-1;
      if(atom >  ibase[try+1]) start = try+1;
   }  while(start <= stop);
   fprintf(out,"Error in GETRES. Couldn't find atom %d\n",atom);
   die();
 
   return(0);
}
 
