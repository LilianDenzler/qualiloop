#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
void SetupSeq(
int nchain,
char *sequence
)
{
int i,len;
 
   /* Get the length of the chain */
   len = strlen(sequence);
 
   /* Copy in the residues */
   strncpy(pstruct.resnme[values.nres++],"NTER",4);
   for(i=0; i<len; i++)
   {
      onethr(sequence[i],pstruct.resnme[values.nres]);
      abmpad(pstruct.resnme[values.nres++],4);
   }
   strncpy(pstruct.resnme[values.nres++],"CTER",4);
 
   /* Update nictot */
   pstruct.segndx[i][0] = values.nres;
   pstruct.segndx[i][1] = 0;
   pstruct.segndx[i][2] = 0;
   pstruct.segndx[i][3] = 0;
   pstruct.segndx[i][4] = 0;
   pstruct.segndx[i][5] = 0;
   pstruct.segndx[i][6] = 0;
   pstruct.segndx[i][7] = 0;
   pstruct.segndx[i][8] = 0;
   pstruct.segndx[i][9] = 0;
}
 
