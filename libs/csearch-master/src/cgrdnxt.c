#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Read a coordinate set from the CG file */
/* com omupd bnj 28/02/92  remove \n from scanf control string... */
 
void cgrdnxt (
FILE *fp,
int *consp,
int ncons,
float *tote,
int *eof
)
{
float x,y,z;
int i;
char InpBuffer[80];

   *eof = 0;
 
   for(i=0;i<ncons;i++)
   {
      /* if(!fscanf(fp,"%f %f %f",&x,&y,&z)) */
      if (!fgets(InpBuffer,80,fp))
      {
         *eof = 1;
         return;
      }
      sscanf(InpBuffer,"%f %f %f",&x,&y,&z);
      coords.xcart[consp[i]-1] = x;
      coords.ycart[consp[i]-1] = y;
      coords.zcart[consp[i]-1] = z;
   }
   /* if(!fscanf(write_fp,"%f",tote)) */
   if (!fgets(InpBuffer,80,fp))
   {
      *eof = 1;
      return;
   }
   sscanf(InpBuffer,"%f",tote);
   return;
}
 
