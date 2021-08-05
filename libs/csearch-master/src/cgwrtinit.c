#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"

/* Write header info to the CG file */
 
void cgwrtinit (
int *consp,
int ncons
)
{
   int i, j;
 
   fprintf(write_fp,"%5d %5d %5d\n",values.nres,values.natoms,ncons);
 
#ifdef SHOWALL
   for(i=0,j=0;i<=values.nres;i++,j++)
   {
      fprintf(write_fp,"%5d ",pstruct.lstatm[i]);
      if(j>11)
      {
         j = -1;
         fprintf(write_fp,"\n");
      }
   }
   if(j) fprintf(write_fp,"\n");
 
   for(i=0,j=0;i<values.nres;i++,j++)
   {
      char temp[8];
      strncpy(temp,pstruct.resnme[i],4);
      temp[4] = '\0';
      fprintf(write_fp,"%4s ",temp);
      if(j>13)
      {
         j = -1;
         fprintf(write_fp,"\n");
      }
   }
   if(j) fprintf(write_fp,"\n");
 
   for(i=0,j=0;i<values.natoms;i++,j++)
   {
      char temp[8];
      strncpy(temp,pstruct.atmnme[i],4);
      temp[4] = '\0';
      fprintf(write_fp,"%4s ",temp);
      if(j>13)
      {
         j = -1;
         fprintf(write_fp,"\n");
      }
   }
   if(j) fprintf(write_fp,"\n");
#endif
 
   for(i=0,j=0;i<ncons;i++,j++)
   {
      fprintf(write_fp,"%5d ",consp[i]);
      if(j>11)
      {
         j = -1;
         fprintf(write_fp,"\n");
      }
   }
   if(j) fprintf(write_fp,"\n");
}
 
