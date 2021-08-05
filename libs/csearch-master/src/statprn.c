#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* ACRM 13.06.91 Rewrite in C
   On the VAX, if running a batch job will convert count and the
   first N integers of a to strings and sets the process name to
   the concatenation of these strings */
 
#ifdef VAX
#include <descrip.h>
#define MAKE_DESCRIP(x,y) x.dsc$w_length=strlen(y),x.dsc$a_pointer=(y)
#endif

void statprn(
int count,
int *a,
int n
)
{
 
/*
 
static int omfirst;
static int ombatch;
char buffer[80];
char *buff;
int  i;
 
   omfirst = TRUE;
 
   $DESCRIPTOR(string_desc," ");
 
   buff = (char *)alloc(16*sizeof(char));
 
   if(omfirst)
   {
      omfirst = FALSE;
      ombatch = batch_job();
   }
   if(ombatch)
   {
      buff[0] = '\0';
      sprintf(buff,"%d ",count);
      for(i = 0; i<n; i++)
      {
         sprintf(buffer,"%d",a[i]);
         buff = strcat(buff,buffer);
         buff = strcat(buff," ");
      }
 
      MAKE_DESCRIP(string_desc,buff);
      if(string_desc.dsc$w_length > 15)
         string_desc.dsc$w_length = 15;
      sys$setprn(&string_desc);
   }
   free(buff);
 
*/
}
