#include "ProtoTypes.h"
#include "CongenProto.h"
 
void setchncls(
float *maxdt,
float *maxg,
float *maxevdw,
logical *cistrans,
int *startres,
logical *ok
)
{
   /* maxdt always has its default */
   *maxdt = cg.maxdt_def;
   /* maxg always has its default */
   *maxg =  100.0;
   /* cistrans is always true */
   *cistrans = f77_true;
 
   sscanf(g_strparam[0],"%f",maxevdw);
 
   padspace(g_strparam[1],4);
   padspace(g_strparam[2],4);
   *startres = get_resnum(g_strparam[1],g_strparam[2]);
 
   if(*startres == 0)
   {
      fprintf(out,"Error in setchncls() -- Starting residue \
incorrectly specified or ommited.\n");
      *ok = f77_false;
   }
}
 
