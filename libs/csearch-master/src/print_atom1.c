#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Converts one atom, iat, into segid res resid iupac and concatenates
     the result onto st. */
 
void print_atom1(
char st[],
int iat
)
{
char buf[21];
int ires,iseg;
 
   if (iat)
   {
      ires = getres((short)iat,pstruct.lstatm,values.nres);
      iseg = getseg(&ires,pstruct.segndx,&values.nsegs);
      sprintf(buf,"%.4s %.4s %.4s %.4s ",
              pstruct.segid[iseg-1],pstruct.resnme[ires-1],
              pstruct.resid[ires-1],pstruct.atmnme[iat-1]);
   }
   else
      strcpy(buf,"       No Atom      ");
   strcat(st,buf);
}
 
