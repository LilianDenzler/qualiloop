#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Read in an atom type and set its mass
   Syntax: TYPE atom-type-code atom-type-name mass

   04.08.92 Some rewriting.   By:   ACRM.
*/
 
void RTFType(
int *NumTypeCodes
)
{
   int AtomTypeCode,
       i;
 
   abmpad(ResTop_strparam[1],4);
 
   sscanf(ResTop_strparam[0],"%d",&AtomTypeCode);

   if(AtomTypeCode < 1 || AtomTypeCode > maxatt)
   {
      fprintf(out,"Error==> Bad atom type code (%d) in residue topology \
ignored\n",AtomTypeCode);
      ResTop_errcnt++;
   }
   else
   {
      *NumTypeCodes = maxi2(AtomTypeCode, *NumTypeCodes);

      /* See if this atom type code has appeared before */
      for(i=0; i<*NumTypeCodes; i++)
      {
         if(!strncmp(restop.acodes[i],ResTop_strparam[1],4))
         {
            fprintf(out,"Error==> Duplicate atom type name in residue topology.\n");
            ResTop_errcnt++;
         }
      }

      strncpy(restop.acodes[AtomTypeCode-1],ResTop_strparam[1],4);
      abmpad(restop.acodes[AtomTypeCode-1],4);

      sscanf(ResTop_strparam[2],"%f",&restop.atmmas[AtomTypeCode-1]);

      if(restop.atmmas[AtomTypeCode-1] <= 0.0)
      {
         fprintf(out,"Error==> Non-positive mass specified in residue topology.\n");
         ResTop_errcnt++;
      }
   }
}
 
