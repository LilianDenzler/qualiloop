#include "ProtoTypes.h"
#include "CongenProto.h"
 
 
/* Clear the atoms required for a CONGEN build */
/* 24.03.94 ACRM Bug fix for PRO at start of loop */
 
void clear_atoms(
char chain[],
char res1[],
char res2[]
)
{
   int atom_start,
       atom_stop,
       atom_ca,
       i,
       count;
 
   if((atom_start = get_atnum(chain,res1,"H   "))==(-1))
   {
      fprintf(out,"Searching for CA...\n");
      atom_start = get_atnum(chain,res1,"CA  ");
   }

   atom_stop  = get_atnum(chain,res2,"C   ");
   atom_ca    = get_atnum(chain,res2,"CA  ");

   if(atom_start==(-1) || atom_stop==(-1) || atom_ca==(-1))
   {
      fprintf(out,"Warning: Atoms not found for clear command\n");
      return;
   }
 
   count = 0;
 
   for(i=atom_start; i<=atom_stop; i++)
      if(i != atom_ca && i != atom_stop)
      {
         coords.xcart[i-1] = coords.ycart[i-1] = coords.zcart[i-1] = 9999.0;
         count++;
      }
 
   fprintf(out,"%d coordinates have been initialised for rebuilding.\n",count);
}
 
