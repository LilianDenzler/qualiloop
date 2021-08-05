#include "ProtoTypes.h"
#include "CongenProto.h"
 
void prtglobopt(
void
)
{
   fprintf(out,"GLYEMAX = %.2g  ALAEMAX = %.2g  PROEMAX = %.2g  \
ERINGPRO = %.2g\n",cg.glyemax,cg.alaemax,cg.proemax,cg.eringpro);
   if (cg.ignore_evdw)
      fprintf(out,"Van der Waals energies will be ignored in energy computations\n");
   if (cg.save_coor)
      fprintf(out,"Coordinates will not be modified.\n");
   else
      fprintf(out,"Coordinates will be set to last generated conformation.\n");
   fprintf(out,"Debug = %d\n",dbg.cgen);
   fprintf(out,"Search will %sbe optimized for FIRST and ALL sidechain\n",
          cg.sidehits_opt ? "" : "NOT ");
   fprintf(out,"options when sidechain contacts occur.\n");
}
 
