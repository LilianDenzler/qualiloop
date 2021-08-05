#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Constructs the atom data structure for an N terminal HT1. */
 
void fillatm_ht1(
struct atom **atp,
int resno,
logical forward
)
{
   int ires,firstat,lastat,fp1,lp1,*antep,iseg,segfirst,seglast;
 
   ires = resno - 1;
   *atp = (struct atom *) alloc(sizeof(struct atom));
   (*atp)->atomno = 0;
   clear_atom(*atp);
   firstat = pstruct.lstatm[ires] + 1;
   lastat = pstruct.lstatm[ires+1];
   fp1 = lastat + 1;
   lp1 = pstruct.lstatm[ires+2];
   iseg = getseg(&resno,pstruct.segndx,&values.nsegs);
   segfirst = pstruct.segndx[iseg-1][0] + 1;
   seglast = pstruct.segndx[iseg][0];
 
   if (forward)
      prdie("Error in fillatm_ht1 -- Forward direction not supported.\n");
 
   (*atp)->atomno = srwdbd(pstruct.atmnme,&firstat,&lastat,"HT1 ");
   if ((*atp)->atomno == 0)
   {
      free(*atp);
      *atp = null;
   }
   else
   {
      antep = &(*atp)->ante[0];
      *antep = srwdbd(pstruct.atmnme,&fp1,&lp1,"C   ");
      *++antep = srwdbd(pstruct.atmnme,&fp1,&lp1,"CA  ");
      *++antep = srwdbd(pstruct.atmnme,&fp1,&lp1,"NT  ");
      fillbndang(*atp);
   }
}
 
