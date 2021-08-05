#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Sets the print flag.
   Syntax: PRINT ON
                 OFF

   02.08.92 Some rewriting.   By:   ACRM.
*/
 
void RTFPrint(
int *PrintFlag
)
{
   if(!strncmp(ResTop_strparam[0],"ON",2))
      *PrintFlag = TRUE;
   else if(!strncmp(ResTop_strparam[0],"OFF",3))
      *PrintFlag = FALSE;
   else
   {
      fprintf(out,"Warning==> Unknown PRINT option in Residue Topology.\n");
      ResTop_wrncnt++;
   }
}
 
