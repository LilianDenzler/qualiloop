#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Copies between entries in the best coordinate arrays. */
 
void cpftbest(
struct sideres *srp,
int fromind,
int toind
)
{
   int *ip;
   float *oxp,*nxp;
   float *oyp,*nyp;
   float *ozp,*nzp;
 
   oxp = srp->bestxpp[fromind];
   oyp = srp->bestypp[fromind];
   ozp = srp->bestzpp[fromind];
   nxp = srp->bestxpp[toind];
   nyp = srp->bestypp[toind];
   nzp = srp->bestzpp[toind];
   for (ip = srp->atomp; *ip; ip++)
   {
      *nxp++ = *oxp++;
      *nyp++ = *oyp++;
      *nzp++ = *ozp++;
   }
}
 
