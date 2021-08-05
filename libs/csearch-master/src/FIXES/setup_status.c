#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Sets up the status degree of freedom. 
   24.05.94 Initialise desc->file_count    By: ACRM
            Corrected desc->nflush to 2
*/
 
void setup_status(
struct status_d *desc,
logical *okp
)
{
   logical ok,changed;
   int maxflush;
 
   desc->setprn = false;
   desc->fileunit = -1;
   desc->nflush = 2;
 
   maxflush = 2;
 
   /* ACRM added casts */
   desc->flushfreq = (int *)alloc(sizeof(int)*maxflush);
   desc->flushunit = (int *)alloc(sizeof(int)*maxflush);
   /* ACRM not sure why this isn't sizeof(int) !!! */
   desc->flush_time = (int *)alloc(sizeof(desc->flush_time)* maxflush);
   desc->filemx = 256;
   desc->filename = alloc(desc->filemx+1);
   strcpy(desc->filename,g_strparam[0]);
   desc->filel = strlen(g_strparam[0]);
   ok = true;
   changed = false;

   /* ACRM 24.05.94 Added this initialisation                           */
   desc->file_count = 0;
 
   /* We don't allow desc->setprn to be true at present */
   /* status_fp is the status file rather than desc->fileunit */
 
   /* Everything will be flushed, etc. every 10 cycles */
   desc->filefreq = 10;
 
   /* We will flush the following units */
   desc->flushunit[0] = (int) out;       /* Cast file pointer to int */
   desc->flushunit[1] = (int) write_fp;  /* Cast file pointer to int */
   /* Every 10 cycles */
   desc->flushfreq[0] = 10;
   desc->flushfreq[1] = 10;
 
   *okp = ok;
}
 
