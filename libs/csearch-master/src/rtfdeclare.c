#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Declare a name-value equivalent
   Syntax: DECLARE name value
   
   02.08.92 Some rewriting.   By:   ACRM
*/
 
void RTFDeclare(
char TableKey[MXTABL][4],
int  *TableReplace,
int  *TableSize
)
{
   int i, bad=FALSE;
 
   abmpad(ResTop_strparam[0],4);
 
   /* Check the symbol doesn't exist already */
   for(i=0; i<*TableSize; i++)
   {
      if(!strncmp(TableKey[i],ResTop_strparam[0],4))
      {
         fprintf(out,"Error==> Symbol redeclared in residue topology\n");
         ResTop_errcnt++;
         bad = TRUE;
         break;
      }
   }
   
   if(!bad)
   {
      if(++(*TableSize) >= MXTABL)
      {
          fprintf(out,"Error==> Too many declarations in residue topology. \
Only %d allowed.\n",MXTABL);
          die();
      }
      strncpy(TableKey[*TableSize-1],ResTop_strparam[0],4);
      sscanf(ResTop_strparam[1],"%d",&TableReplace[*TableSize-1]);
   }
}
 
