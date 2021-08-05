#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Read a PIR file containing multiple chains of up to maxaa amino acids.
   Each chain is returned in seqs[].
   The number of chains is returned by the routine.
 
   In:      FILE *fp       File pointer
            int  maxaa     Max number of aa's in a chain
   Out:     char **seqs    Array of pointers to sequences
   Return:  int            Number of chains.
*/
 
int ReadPIR(
FILE *fp,
int  maxaa,
char *seqs[]
)
{
   char *buffer;
   int  rescount = 0,
        chain = 0;
 
   /* Allocate space for the sequence */
   buffer = (char *)malloc((maxaa+1) * sizeof(char));
   if(!buffer) prdie("ReadPIR(): Unable to allocate buffer space.");
 
   /* Read header lines from the file */
   fgets(buffer,maxaa,fp);
   fgets(buffer,maxaa,fp);
 
   /* Now loop through to get the sequence */
   while(rescount<maxaa && !feof(fp))
   {
      int ch;
 
      /* Get a character */
      ch = getc(fp);
      if(ch==EOF) break;
      buffer[rescount] = ch;
 
      if(isalpha(buffer[rescount]))
      {
         /* If it's an alpha character, then toupper() it and
            increment the counter.
         */
         buffer[rescount] = toupper(buffer[rescount]);
         rescount++;
      }
      else if(buffer[rescount] == '*')
      {
         /* If it's a star, then it's the end of a chain, so copy the chain */
         buffer[rescount] = '\0';
         seqs[chain] = (char *)malloc((rescount+2)*sizeof(char));
         if(!seqs[chain]) prdie("ReadPIR(): Unable to allocate sequence space.");
         strcpy(seqs[chain],buffer);
         chain++;
         rescount=0;
      }
   }
 
   /* Check to see if the last chain ended without a * */
   if(rescount)
   {
      buffer[rescount] = '\0';
      seqs[chain] = (char *)malloc((rescount+2)*sizeof(char));
      if(!seqs[chain]) prdie("ReadPIR(): Unable to allocate sequence space.");
      strcpy(seqs[chain],buffer);
      chain++;
   }
   free(buffer);
 
   return(chain);
}
 
