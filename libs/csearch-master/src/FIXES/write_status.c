#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* One day this could be improved to keep track of a version number
   so that two status files are kept at a time 

   24.05.94 Now reopens the file if it was already open and closed.
            Uses the correct filename (!)
            Added () round write_count
            Does a rename() to keep the last version
            By: ACRM
*/
 
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
   char day[20],
        file[160],
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

   /* ACRM 24.05.94 Added:                                              */
   sprintf(file,"%s.old",fst);
   rename(fst,file);

   (*write_count)++;

/* ACRM 24.05.94:
// if((status_fp = fopen(file,"w"))==NULL)
*/
   if((status_fp = fopen(fst,"w"))==NULL)
   {
      fprintf(out,"Warning: Unable to open status file\n");
      return;
   }

   /* get the date and time */
   GetDate(day);
   GetTime(time);
 
   /* Write to the status file and close it */
   fprintf(status_fp,"Time: %s %s.  Leafnum is %d\nIters = ",
           day,time,count);
   for(i=0;i<niter;i++)
      fprintf(status_fp,"%d ",iters[i]);
   fprintf(status_fp,"\n");
   fclose(status_fp);
   status_fp = NULL;

   return;
}
 
