#include "ProtoTypes.h"
#include "CongenProto.h"

/***********************************************************************/ 
/*>void clear_sc(char chain[],char res1[],char res2[])
   ---------------------------------------------------
   Clear sidechain atoms in a specified residue range. This is the
   *correct* way to do it. Replaces OML's clear_side() routine which
   doesn't handle chain labels properly.

   28.03.94 Original    By: ACRM 
*/
void clear_sc(char chain[],char res1[],char res2[])
{
   int atom_start,
       atom_stop,
       iatom,
       count;
 
   atom_start = get_atnum(chain,res1,"N   ");
   atom_stop  = get_atnum(chain,res2,"C   ");
 
   if(atom_start == (-1) || atom_stop == (-1))
   {
      fprintf(out,"Error: Residue specified for SCCLEAR not found.\n");
      return;
   }

   count = 0;
 
   for(iatom=atom_start; iatom<=atom_stop; iatom++)
   {
      
      if(strncmp(pstruct.atmnme[iatom-1],"N   ",4) &&
         strncmp(pstruct.atmnme[iatom-1],"H   ",4) &&
         strncmp(pstruct.atmnme[iatom-1],"CA  ",4) &&
         strncmp(pstruct.atmnme[iatom-1],"CB  ",4) &&
         strncmp(pstruct.atmnme[iatom-1],"C   ",4) &&
         strncmp(pstruct.atmnme[iatom-1],"O   ",4))
      {
         coords.xcart[iatom-1] = 
         coords.ycart[iatom-1] = 
         coords.zcart[iatom-1] = (float)9999.0;

         count++;
      }
   }
 
   fprintf(out,"%d coordinates have been initialised for rebuilding.\n",
           count);
}
 
