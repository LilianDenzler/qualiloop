#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* One day this could be improved to keep track of a version number
   so that two status files are kept at a time */
 
void write_status(
int count,
int *iters,
int niter,
int unit,
int *write_count,
char fst[],
int fstmax,
int fstlen
)
{
   char file[256],
        day[20],
        time[20];
   int  exists,
        opened,
        i;
 
   /* If we already have a status file for some reason */
   if(status_fp)
   {
      fclose(status_fp);
      status_fp = NULL;
   }
   else
   {
      *write_count++;
      if((status_fp = fopen(file,"w"))==NULL)
      {
         fprintf(out,"Warning: Unable to open status file\n");
         return;
      }
   }
   /* get the date and time */
   GetDate(day);
   GetTime(time);
 
   /* Write to the status file and close it */
   fprintf(status_fp,"Time: %s %s.  Leafnum is %d\nIters = ",day,time,count);
   for(i=0;i<niter;i++)
      fprintf(status_fp,"%d ",iters[i]);
   fprintf(status_fp,"\n");
   fclose(status_fp);
   status_fp = NULL;
 
   return;
}
 
