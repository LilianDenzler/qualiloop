#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Finds the atom with iupac name in residue ires. It is assumed that iupac
*   is a four character name.
*/
 
int typeinres(
int ires,
char iupac[]
)
{
   int firstat,lastat,i;
 
   firstat = pstruct.lstatm[ires-1]+1;
   lastat  = pstruct.lstatm[ires];
 
   i = srwdbd(pstruct.atmnme,&firstat,&lastat,iupac);
 
   return i;
}
 
