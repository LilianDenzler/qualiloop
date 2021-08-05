#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Constructs the atom data structure for a backbone H. */
 
void fillatm_h(
struct atom **atp,
int resno,
logical forward
)
{
   int ires,firstat,lastat,fm1,lm1,fp1,lp1,*antep,iseg,segfirst,seglast;
 
   ires = resno - 1;
   *atp = (struct atom *) alloc(sizeof(struct atom));
   (*atp)->atomno = 0;
   clear_atom(*atp);
   fm1 = pstruct.lstatm[ires-1]+1;
   lm1 = pstruct.lstatm[ires];
   firstat = lm1 + 1;
   lastat = pstruct.lstatm[ires+1];
   fp1 = lastat + 1;
   lp1 = pstruct.lstatm[ires+2];
   iseg = getseg(&resno,pstruct.segndx,&values.nsegs);
   segfirst = pstruct.segndx[iseg-1][0] + 1;
   seglast = pstruct.segndx[iseg][0];
 
   (*atp)->atomno = srwdbd(pstruct.atmnme,&firstat,&lastat,
                           resno == segfirst+1 ? "HT3 " : "H   ");
   if ((*atp)->atomno == 0)
   {
      free(*atp);
      *atp = null;
   }
   else
   {
      antep = &(*atp)->ante[0];
      if (forward)
      {
         if (!(*antep = srwdbd(pstruct.atmnme,&fm1,&lm1,"CA  ")) &&
                resno == segfirst+1)
            *antep = srwdbd(pstruct.atmnme,&fm1,&lm1,"HT2 ");
         if (!(*++antep = srwdbd(pstruct.atmnme,&fm1,&lm1,"C   ")) &&
                resno == segfirst+1)
            *antep = srwdbd(pstruct.atmnme,&fm1,&lm1,"HT1 ");
         *++antep = srwdbd(pstruct.atmnme,&firstat,&lastat,
                           resno==segfirst+1 ? "NT  " : "N   ");
      }
      else
      {
         *antep = srwdbd(pstruct.atmnme,&firstat,&lastat,"C   ");
         *++antep = srwdbd(pstruct.atmnme,&firstat,&lastat,"CA  ");
         *++antep = srwdbd(pstruct.atmnme,&firstat,&lastat,
                           resno==segfirst+1 ? "NT  " : "N   ");
      }
      fillbndang(*atp);
   }
}
 
