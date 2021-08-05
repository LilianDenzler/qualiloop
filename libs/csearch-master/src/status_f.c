#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Status degree of freedom. This routine will set the process name
*   to the concatenation of the leaf number and iteration counts of the first
*   degrees of freedom so that the user can tell how far CGEN has proceeded,
*   write this same information repeatedly to a file, and periodically close
*   and reopen the output file (unit 6).
*
*   24.05.94 Commented properly and removed check on desc->fileunit
*            By: ACRM
*/
 
void status_f(
struct dof *dofp,
struct status_d *desc
)
{
   int curtime;
   int *iters;
   struct dof *p;
   int ni,*ip,i;
 
   /* Increment the call count and get the time                         */
   desc->callcount++;
   curtime = time(null);
   if (dbg.cgen)
      fprintf(out,"status_f: leafnum = %d  time = %d  callcount = %d\n",
             leafnum,curtime,desc->callcount);

   /* Count the number of degrees of freedom (excluding status) and 
      allocate memory for an array of this size
   */
   ni = 0;
   for (p=dof_head; p!=dofp; p=p->next)
      if (p->dof_type != status_t) ni++;
   ip = iters = (int *) alloc(ni*sizeof(int));

   /* Walk through the degrees of freedom and store the iteration counts
      into this new array
   */
   for (p=dof_head; p!=dofp; p=p->next)
      if (p->dof_type != status_t) *ip++ = p->iter;

   /* Set the process name if required                                  */
   if (desc->setprn) statprn(leafnum,iters,ni);

   /* If we have a file unit for this degree of freedom and a frequency
      >= 0 then write status
   */
/* ACRM 24.05.94: 
// if (desc->fileunit >= 0 && desc->filefreq >= 0)
*/
   if (desc->filefreq >= 0)
   {
      if (desc->callcount % desc->filefreq == 0 ||
          curtime-desc->file_time >= 3600)
      {
/* ACRM 24.05.94: Added () round file_count
*/
         write_status(leafnum,iters,ni,desc->fileunit,
                      &(desc->file_count),desc->filename,desc->filemx,
                      desc->filel);
         desc->file_time = curtime;
      }
   }
   
   /* Flush the output files                                            */
   for (i=0; i<desc->nflush; i++)
   {
      if (desc->callcount % desc->flushfreq[i] == 0 ||
         curtime-desc->flush_time[i] >= 3600)
      {
/*       Replaced with C flush
//       flush_output(&desc->flushunit[i]);
*/
         fflush((FILE *)desc->flushunit[i]);
         desc->flush_time[i] = curtime;
      }
   }

   /* Free the array of iteration counts                                */
   free(iters);

   /* Recursive dispatch of next DoF                                    */
   dispatch(dofp->next);
}
 
