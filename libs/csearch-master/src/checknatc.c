#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Check the maximum number of atom types has not been exceeded
   02.08.92 Some recoding.    By:   ACRM
*/
 
void CheckNATC(int *NumAtTypeCodes)
{
   int i;
   
   int LocAtCodes = 0;

   for(i=0; i<*NumAtTypeCodes; i++)
      if(strncmp(restop.acodes[i],BLANK,4)) LocAtCodes++;

   if(LocAtCodes > maxatu)
   {
      fprintf(out,"Error==> Number of atom types (%d) in residue topology \
exceeds maximum (%d)\n",LocAtCodes,maxatu);
   }
}
 
