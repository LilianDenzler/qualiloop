#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Write a coordinate set to the CG file */
 
void cgwrtnxt (
int *consp,
int ncons,
float tote
)
{
   int i;
 
   for(i=0;i<ncons;i++)
   {
      fprintf(write_fp,"%8.3f %8.3f %8.3f\n",
      coords.xcart[consp[i]-1],coords.ycart[consp[i]-1],
      coords.zcart[consp[i]-1]);
   }
   fprintf(write_fp,"%15.5f\n",tote);
}
 
