#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Someday soon, this should freeup storage */
 
void termcgen(
void
)
{
   fprintf(out,"\nA total of %d conformations has been written.\n",confnum);
   if (confnum != leafnum)
      fprintf(out,"A total of %d bottom leaves were hit.\n",leafnum);
   if (cg.save_coor)
   {
/* Following is commented out until HBOND stops generating hydrogen positions.
      int i,i1;
 
      for (i=0; i<cg.nconsp; i++)
      {
         i1 = cg.consp[i] - 1;
         coords.xcart[i1] = cg.savex[i];
         coords.ycart[i1] = cg.savey[i];
         coords.zcart[i1] = cg.savez[i];
      }
*/
      CopyFArray(cg.savex,coords.xcart,values.natoms);
      CopyFArray(cg.savey,coords.ycart,values.natoms);
      CopyFArray(cg.savez,coords.zcart,values.natoms);
      free(cg.savex);
      free(cg.savey);
      free(cg.savez);
   }
#ifdef EVALUATE
   for (p=dof_head; p; p=p->next)
   {
      if (p->dof_type == evaluate_t)
      {
         desc = p->desc;
         if (desc->evalcount > 0)
         {
            fprintf(out,"Minimum evaluation was %g\n",desc->mineval);
            fprintf(out,"Maximum evaluation was %g\n",desc->maxeval);
         }
      }
   }
#endif /* EVALUATE */
}
 
