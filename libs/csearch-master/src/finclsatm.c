#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/*
*   This routine actually calculates the distance from the current
*   residue to the closing atom and adjusts it as specified by the DELTA
*   option in the BACKBONE degree of freedom.
*
*   The routine first determines the list of atoms that must be
*   constructed in order to go from the last atom constructed by this
*   degree of freedom (pointed to by desc) to the closing atom. Then,
*   the routine finds the maximum distance that can be spanned by these
*   atoms. The adjustments to the bond angles necessitated by the chain
*   closure procedure is handled by constructing a list of bond angles
*   that can be stretched and matching against it when the model is
*   constructed.
*/
 
void finclsatm(
struct backbone_d *desc
)
{
   int natom;
   int *catom;
   int **cangle;
   int nclschn,ncangle,iangle;
   logical stretch_all;
   float d,*cmaxdt;
 
   if (desc->closing_atom > 0)
   {
      if (desc->qdelta || desc->closing_distance == 0.0)
      {
         catom = (int *) alloc(values.natoms * sizeof(int));
         stretch_all = desc->maxdt != -1.0;
         if (stretch_all)
         {
            cmaxdt = (float *) alloc(sizeof(float));
            *cmaxdt = desc->maxdt;
         }
         else
         {
            nclschn = 0;
            walk_dofs(if (dofp->dof_type == chain_closure_t) nclschn++; );
            if (nclschn == 0)
            {
               stretch_all = true;
               cmaxdt = (float *) alloc(sizeof(float));
               *cmaxdt = cg.maxdt_def * dtorad;
               fprintf(out,"\nNote from finclsatm -- \
No chain closures found\n");
               fprintf(out,"MAXDT for CLSA for residue %.4s %.4s will be \
set to %g.\n", pstruct.resnme[desc->resno-1], pstruct.resid[desc->resno-1],
                      (*cmaxdt)/dtorad);
            }
            else
            {
               int i;

               ncangle = nclschn * 9;
               cmaxdt = (float *) alloc(ncangle * sizeof(float));
               cangle = (int **)alloc(ncangle * sizeof(int *));
               for(i=0; i<ncangle; i++)
               {
                  cangle[i] = (int *)alloc(3 * sizeof(int));
               }
               iangle = 0;
               walk_dofs(
                  if (dofp->dof_type == chain_closure_t)
                  {
                     /* ACRM added cast */
                     adccangs((int **)cangle,cmaxdt,&iangle,
                        (struct chain_closure_d *)dofp->desc);
                  }
               );
               if (iangle != ncangle)
               {
                  fprintf(out,"\nError in finclsatm -- \
iangle+1 != ncangle\n");
                  die();
               }
            }
         }
         if (desc->forward)
            cnsfblist(desc,catom,&natom);
         else
            cnsbblist(desc,catom,&natom);
 
         if (natom == 0)
            desc->closing_distance = 0.0;
         else
         {
            d = desc->closing_distance;
            cclsngdst(desc,catom,natom,stretch_all,
                                       cangle,cmaxdt,ncangle);
            if (desc->qdelta) desc->closing_distance += d;
         }
         free(catom);
         if (!stretch_all) free(cangle);
         free(cmaxdt);
      }
      fprintf(out,"For backbone degree of freedom %.4s %.4s\n",
             pstruct.resnme[desc->resno-1],pstruct.resid[desc->resno-1]);
      fprintf(out,"the final closing distance is %g\n",
             desc->closing_distance);
   }
}
 
