#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Returns the RTF TOPOLOGY type depending on the first letter of the
   keyword.
                    { PROT }
                    { HPRO }
   Syntax: TOPOLOGY { ALLH }
                    { DNA  }
                    { UNKN }
   02.08.92 Some rewriting.   By:   ACRM.
*/
 
int RTFTopology(
void
)
{
   int RTFType;
   
   switch(ResTop_strparam[0][0])
   {
   case 'P':
   case 'p':
      RTFType = 1;
      break;
   case 'H':
   case 'h':
      RTFType = 2;
      break;
   case 'A':
   case 'a':
      RTFType = 3;
      break;
   case 'D':
   case 'd':
      RTFType = 4;
      break;
   case 'U':
   case 'u':
      RTFType = 0;
      break;
   default:
      fprintf(out,"Warning==> Unrecognised Residue topology type.\n");
      fprintf(out,"           UNKNOWN will be used.\n");
      ResTop_wrncnt++;
      RTFType = 0;
      break;
   }
   
   return(RTFType);
}
 
